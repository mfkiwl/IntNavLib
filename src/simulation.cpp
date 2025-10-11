#include "simulation.h"

namespace intnavlib {

ImuMeasurements kinematicsEcef(const NavSolutionEcef & new_nav,
                                const NavSolutionEcef & old_nav) {
    
    double tor_i = new_nav.time - old_nav.time;

    // Init measurements to 0
    ImuMeasurements true_imu_meas;
    true_imu_meas.time = new_nav.time;
    true_imu_meas.f = Eigen::Vector3d::Zero();
    true_imu_meas.omega = Eigen::Vector3d::Zero();

    if (tor_i > 0.0) {

    // From (2.145) determine the Earth rotation over the update interval
    // C_Earth = C_e_i' * old_C_e_i

    double alpha_ie = kOmega_ie * tor_i;

    double cos_alpha_ie = cos(alpha_ie);
    double sin_alpha_ie = sin(alpha_ie);

    Eigen::Matrix3d C_Earth;
    C_Earth << cos_alpha_ie, sin_alpha_ie, 0.0,
                -sin_alpha_ie, cos_alpha_ie, 0.0,
                    0.0,             0.0,  1.0;

    // Obtain coordinate transformation matrix from the old attitude (w.r.t.
    // an inertial frame) to the new (compensate for earth rotation)
    Eigen::Matrix3d C_old_new = new_nav.C_b_e.transpose() * C_Earth * old_nav.C_b_e;

    // Calculate the approximate angular rate w.r.t. an inertial frame
    Eigen::Vector3d alpha_ib_b;
    alpha_ib_b(0) = 0.5 * (C_old_new(1,2) - C_old_new(2,1));
    alpha_ib_b(1) = 0.5 * (C_old_new(2,0) - C_old_new(0,2));
    alpha_ib_b(2) = 0.5 * (C_old_new(0,1) - C_old_new(1,0));

    // Calculate and apply the scaling factor
    double temp = acos(0.5 * (C_old_new(0,0) + C_old_new(1,1) + C_old_new(2,2) - 1.0));
    if (temp > 2.0e-5) // scaling is 1 if temp is less than this
        alpha_ib_b = alpha_ib_b * temp/sin(temp);
    
    // Calculate the angular rate
    true_imu_meas.omega = alpha_ib_b / tor_i;

    // Calculate the specific force resolved about ECEF-frame axes
    // From (5.36)
    Eigen::Vector3d kOmega_ie_vec;
    kOmega_ie_vec << 0.0 , 0.0 , kOmega_ie;
    Eigen::Vector3d f_ib_e = ((new_nav.v_eb_e - old_nav.v_eb_e) / tor_i) - gravityEcef(old_nav.r_eb_e)
        + 2.0 * skewSymmetric(kOmega_ie_vec) * old_nav.v_eb_e;

    // Calculate the average body-to-ECEF-frame coordinate transformation
    // matrix over the update interval using (5.84) and (5.85)

    double mag_alpha = alpha_ib_b.norm();
    Eigen::Matrix3d Alpha_ib_b = skewSymmetric(alpha_ib_b);

    Eigen::Vector3d alpha_ie_vec;
    alpha_ie_vec << 0.0 , 0.0 , alpha_ie;   

    Eigen::Matrix3d ave_C_b_e;
    if (mag_alpha>1.0e-8) {
        ave_C_b_e = old_nav.C_b_e * 
            (Eigen::Matrix3d::Identity() + 

            ((1.0 - cos(mag_alpha)) / pow(mag_alpha,2.0)) *
            Alpha_ib_b + 
            
            ((1.0 - sin(mag_alpha) / mag_alpha) / pow(mag_alpha,2.0)) * 
            Alpha_ib_b * Alpha_ib_b) - 

            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    else { // Approximate if angle small enough (sum not multiply)
        ave_C_b_e = old_nav.C_b_e -
            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    
    // Transform specific force to body-frame resolving axes using (5.81)
    // plain inverse as matrix is 3x3
    // true_imu_meas.f = ave_C_b_e.inverse() * f_ib_e;
    true_imu_meas.f = ave_C_b_e.colPivHouseholderQr().solve(f_ib_e);

    }

    return true_imu_meas;

}


ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas,
                        const ImuMeasurements & old_imu_meas,
                        const ImuErrors & imu_errors,
                        const double & tor_i,
                        std::mt19937 & gen) {
    
    ImuMeasurements imu_measurements;

    imu_measurements.time = true_imu_meas.time;

    // Init noises
    Eigen::Vector3d accel_noise = Eigen::Vector3d::Zero();
    Eigen::Vector3d gyro_noise = Eigen::Vector3d::Zero();

    // Normal distribution init
    std::normal_distribution<double> randn(0.0, 1.0); 

    if(tor_i > 0.0) {
        accel_noise << randn(gen) * imu_errors.accel_noise_root_psd / sqrt(tor_i), 
                        randn(gen) * imu_errors.accel_noise_root_psd / sqrt(tor_i),
                        randn(gen) * imu_errors.accel_noise_root_psd / sqrt(tor_i);  
        gyro_noise << randn(gen) * imu_errors.gyro_noise_root_psd / sqrt(tor_i),
                        randn(gen) * imu_errors.gyro_noise_root_psd / sqrt(tor_i),
                        randn(gen) * imu_errors.gyro_noise_root_psd / sqrt(tor_i);
    }

    // Calculate accelerometer and gyro outputs using (4.16) and (4.17)

    // Specific force
    Eigen::Vector3d uq_f_ib_b = imu_errors.b_a + (Eigen::Matrix3d::Identity() + imu_errors.M_a) * true_imu_meas.f 
                                + accel_noise;

    // Angular velocity
    Eigen::Vector3d uq_omega_ib_b = imu_errors.b_g + (Eigen::Matrix3d::Identity() + imu_errors.M_g) * true_imu_meas.omega 
                                    + imu_errors.G_g * true_imu_meas.f 
                                    + gyro_noise;

    // Quantize accelerometer outputs
    if (imu_errors.accel_quant_level > kEpsilon) {
        Eigen::Vector3d temp = (uq_f_ib_b + old_imu_meas.quant_residuals_f) / imu_errors.accel_quant_level;
        for(size_t i = 0; i < 3; i++) temp(i) = round(temp(i));
        imu_measurements.f = imu_errors.accel_quant_level * temp;
        imu_measurements.quant_residuals_f = uq_f_ib_b + old_imu_meas.quant_residuals_f -
            imu_measurements.f;
    }
    else{
        imu_measurements.f = uq_f_ib_b;
        imu_measurements.quant_residuals_f = Eigen::Vector3d::Zero();
    }

    // Quantize gyro outputs
    if (imu_errors.gyro_quant_level > kEpsilon) {
        Eigen::Vector3d temp = (uq_omega_ib_b + old_imu_meas.quant_residuals_omega) / imu_errors.gyro_quant_level;
        for(size_t i = 0; i < 3; i++) temp(i) = round(temp(i)); 
        imu_measurements.omega = imu_errors.gyro_quant_level * temp;
        imu_measurements.quant_residuals_omega = uq_omega_ib_b + old_imu_meas.quant_residuals_omega -
            imu_measurements.omega;
    }
    else {
        imu_measurements.omega = uq_omega_ib_b;
        imu_measurements.quant_residuals_omega = Eigen::Vector3d::Zero();
    }

    // // No quantization
    // imu_measurements.f = uq_f_ib_b;
    // imu_measurements.omega = uq_omega_ib_b;

    return imu_measurements;
}

PosMeasEcef genericPosSensModel(const NavSolutionEcef & true_nav, 
                                const double & pos_sigma,
                                std::mt19937 & gen){

    PosMeasEcef pos_meas;
    pos_meas.time = true_nav.time;

    // nomally distributed error
    std::normal_distribution<double> randn(0.0, pos_sigma);

    Eigen::Vector3d pos_error;
    pos_error << randn(gen), randn(gen), randn(gen);

    Eigen::Matrix3d cov_mat = Eigen::Matrix3d::Identity() * pos_sigma * pos_sigma;

    // add error
    pos_meas.r_eb_e = true_nav.r_eb_e + pos_error;
    pos_meas.cov_mat = cov_mat;

    return pos_meas;
}

PosRotMeasEcef genericPosRotSensModel(const NavSolutionEcef & true_nav, 
                                        const double & pos_sigma,
                                        const double & rot_sigma,
                                        std::mt19937 & gen){
    PosRotMeasEcef pos_rot_meas;
    pos_rot_meas.time = true_nav.time;

    // nomally distributed error
    std::normal_distribution<double> randn_pos(0.0, pos_sigma);
    std::normal_distribution<double> randn_rot(0.0, rot_sigma);

    Eigen::Vector3d pos_error;
    pos_error << randn_pos(gen), randn_pos(gen), randn_pos(gen);

    Eigen::Vector3d rot_error;
    rot_error << randn_rot(gen), randn_rot(gen), randn_rot(gen);
    Eigen::Matrix3d C_b_b = eulerToDcm(rot_error); // perturbation Rot mat

    Eigen::Matrix<double,6,6> cov_mat = Eigen::Matrix<double,6,6>::Identity();
    cov_mat.block<3,3>(0,0) *= pos_sigma * pos_sigma;
    cov_mat.block<3,3>(3,3) *= rot_sigma * rot_sigma;

    // add error
    pos_rot_meas.r_eb_e = true_nav.r_eb_e + pos_error;
    pos_rot_meas.C_b_e = true_nav.C_b_e * C_b_b;
    pos_rot_meas.cov_mat = cov_mat;

    return pos_rot_meas;
}

SatPosVel satellitePositionsAndVelocities(const double & time, 
                                        const GnssConfig& gnss_config) {

    SatPosVel gnssPosVel;

    // Convert inclination angle to radians
    double inclination = gnss_config.inclination * kDegToRad; 

    // Determine orbital angular rate
    double omega_is = std::sqrt(kGravConst / std::pow(gnss_config.r_os, 3));

    // Determine constellation time
    double const_time = time + gnss_config.const_delta_t;

    // Resize output matrices
    int no_sat = static_cast<int>(gnss_config.no_sat);

    gnssPosVel.sat_r_es_e = Eigen::MatrixXd(no_sat, 3);
    gnssPosVel.sat_v_es_e = Eigen::MatrixXd(no_sat, 3);

    // Loop over satellites
    for (int j = 0; j < no_sat; ++j) {
        // Argument of latitude
        double u_os_o = 2 * M_PI * j / no_sat + omega_is * const_time;
        
        // Satellite position in the orbital frame
        Eigen::Vector3d r_os_o;
        r_os_o << gnss_config.r_os * std::cos(u_os_o),
                    gnss_config.r_os * std::sin(u_os_o), 
                    0;

        // Longitude of the ascending node
        double Omega = (M_PI * (j % 6) / 3 + gnss_config.const_delta_lambda * kDegToRad) - kOmega_ie * const_time;

        // ECEF Satellite Position
        gnssPosVel.sat_r_es_e(j, 0) = r_os_o(0) * std::cos(Omega) - r_os_o(1) * std::cos(inclination) * std::sin(Omega);
        gnssPosVel.sat_r_es_e(j, 1) = r_os_o(0) * std::sin(Omega) + r_os_o(1) * std::cos(inclination) * std::cos(Omega);
        gnssPosVel.sat_r_es_e(j, 2) = r_os_o(1) * std::sin(inclination);

        // Satellite velocity in the orbital frame
        Eigen::Vector3d v_os_o;
        v_os_o << -gnss_config.r_os * omega_is * std::sin(u_os_o),
                    gnss_config.r_os * omega_is * std::cos(u_os_o), 0;

        // ECEF Satellite velocity
        gnssPosVel.sat_v_es_e(j, 0) = v_os_o(0) * std::cos(Omega) - v_os_o(1) * std::cos(inclination) * std::sin(Omega) + kOmega_ie * gnssPosVel.sat_r_es_e(j, 1);
        gnssPosVel.sat_v_es_e(j, 1) = v_os_o(0) * std::sin(Omega) + v_os_o(1) * std::cos(inclination) * std::cos(Omega) - kOmega_ie * gnssPosVel.sat_r_es_e(j, 0);
        gnssPosVel.sat_v_es_e(j, 2) = v_os_o(1) * std::sin(inclination);
    }

    return gnssPosVel;
}

GnssMeasurements generateGnssMeasurements(const double & time,
                                        const SatPosVel & gnss_pos_vel,
                                        const NavSolutionNed& true_nav_ned,
                                        const NavSolutionEcef& true_nav_ecef,
                                        const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxGnssSatellites, 1> & gnss_biases, 
                                        const GnssConfig& gnss_config,
                                        std::mt19937 & gen) {

    // nomally distributed error
    std::normal_distribution<double> randn(0.0, 1.0);

    GnssMeasurements gnss_measurements;

    // Initialize number of GNSS measurements
    gnss_measurements.no_meas = 0;

    // Calculate ECEF to NED coordinate transformation matrix
    double cos_lat = std::cos(true_nav_ned.latitude);
    double sin_lat = std::sin(true_nav_ned.latitude);
    double cos_long = std::cos(true_nav_ned.longitude);
    double sin_long = std::sin(true_nav_ned.longitude);
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
             -sin_long,            cos_long,        0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
     
    // Skew symmetric matrix of Earth rate
    Eigen::Matrix3d Omega_ie = skewSymmetric(Eigen::Vector3d(0, 0, kOmega_ie));
       
    // Resize the GNSS measurements matrix to accommodate the maximum possible measurements
    gnss_measurements.meas = Eigen::MatrixXd(gnss_config.no_sat, 8);

    // Loop over satellites
    for (int j = 0; j < gnss_config.no_sat; ++j) {
        // Determine ECEF line-of-sight vector
        Eigen::Vector3d delta_r = gnss_pos_vel.sat_r_es_e.row(j).transpose() - true_nav_ecef.r_eb_e;
        double approx_range = delta_r.norm();
        Eigen::Vector3d u_as_e = delta_r / approx_range;
    
        // Convert line-of-sight vector to NED and determine elevation
        double elevation = -std::asin(C_e_n.row(2).dot(u_as_e));
    
        // Determine if satellite is above the masking angle
        if (elevation >= gnss_config.mask_angle * kDegToRad) {
            // Increment number of measurements
            gnss_measurements.no_meas++;
    
            // Calculate frame rotation during signal transit time
            Eigen::Matrix3d C_e_I;
            C_e_I << 1, kOmega_ie * approx_range / c, 0,
                     -kOmega_ie * approx_range / c, 1, 0,
                     0, 0, 1;

            // Calculate range
            delta_r = C_e_I * gnss_pos_vel.sat_r_es_e.row(j).transpose() - true_nav_ecef.r_eb_e;
            double range = delta_r.norm();
        
            // Calculate range rate
            double range_rate = u_as_e.dot(C_e_I * (gnss_pos_vel.sat_v_es_e.row(j).transpose() + Omega_ie * gnss_pos_vel.sat_r_es_e.row(j).transpose()) - (true_nav_ecef.v_eb_e + Omega_ie * true_nav_ecef.r_eb_e));
    
            // Calculate pseudo-range measurement
            gnss_measurements.meas(gnss_measurements.no_meas-1, 0) = range + gnss_biases(j) + gnss_config.rx_clock_offset + gnss_config.rx_clock_drift * time + gnss_config.code_track_err_sd * randn(gen);
    
            // Calculate pseudo-range rate measurement
            gnss_measurements.meas(gnss_measurements.no_meas-1, 1) = range_rate + gnss_config.rx_clock_drift + gnss_config.rate_track_err_sd * randn(gen);
    
            // Append satellite position and velocity to output data
            gnss_measurements.meas.block<1, 3>(gnss_measurements.no_meas-1, 2) = gnss_pos_vel.sat_r_es_e.row(j);
            gnss_measurements.meas.block<1, 3>(gnss_measurements.no_meas-1, 5) = gnss_pos_vel.sat_v_es_e.row(j);
        }
    }

    // Resize the GNSS measurements matrix to the actual number of measurements
    // Which is <= than n of satellites
    gnss_measurements.meas.conservativeResize(gnss_measurements.no_meas, 8);

    // 6. Set-up measurement noise covariance matrix assuming all measurements
    // are independent and have equal variance for a given measurement type.
    gnss_measurements.cov_mat =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites>::Identity(2*gnss_measurements.no_meas, 2*gnss_measurements.no_meas);

    // Ranges
    gnss_measurements.cov_mat.block(0,0,gnss_measurements.no_meas,gnss_measurements.no_meas) 
                *= pow(gnss_config.pseudo_range_sd,2.0);
    // Range rates
    gnss_measurements.cov_mat.block(gnss_measurements.no_meas,gnss_measurements.no_meas,gnss_measurements.no_meas,gnss_measurements.no_meas)
                *= pow(gnss_config.range_rate_sd,2.0);

    return gnss_measurements;
}

Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxGnssSatellites, 1>
initializeGnssBiases(const NavSolutionEcef & true_nav_ecef,
                                    const NavSolutionNed & true_nav_ned,
                                    const SatPosVel & gnss_pos_vel,
                                    const GnssConfig& gnss_config,
                                    std::mt19937 & gen) {

    // nomally distributed error
    std::normal_distribution<double> randn(0.0, 1.0);

    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxGnssSatellites, 1> gnss_biases(gnss_config.no_sat);

    // Calculate ECEF to NED coordinate transformation matrix
    double cos_lat = std::cos(true_nav_ned.latitude);
    double sin_lat = std::sin(true_nav_ned.latitude);
    double cos_long = std::cos(true_nav_ned.longitude);
    double sin_long = std::sin(true_nav_ned.longitude);
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
             -sin_long,            cos_long,        0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;

    // Loop over satellites
    for (int j = 0; j < gnss_config.no_sat; ++j) {
        // Determine ECEF line-of-sight vector
        Eigen::Vector3d delta_r = gnss_pos_vel.sat_r_es_e.row(j).transpose() - true_nav_ecef.r_eb_e;
        Eigen::Vector3d u_as_e = delta_r / delta_r.norm();
    
        // Convert line-of-sight vector to NED and determine elevation
        double elevation = -std::asin(C_e_n.row(2).dot(u_as_e));
    
        // Limit the minimum elevation angle to the masking angle
        elevation = std::max(elevation, gnss_config.mask_angle * kDegToRad);
    
        // Calculate ionosphere and troposphere error standard deviations
        double iono_SD = gnss_config.zenith_iono_err_sd / std::sqrt(1 - 0.899 * std::cos(elevation) * std::cos(elevation));
        double trop_SD = gnss_config.zenith_trop_err_sd / std::sqrt(1 - 0.998 * std::cos(elevation) * std::cos(elevation));
    
        // Determine range bias
        gnss_biases(j) = gnss_config.sis_err_sd * randn(gen) + iono_SD * randn(gen) + trop_SD * randn(gen);
    }

    return gnss_biases;
}


};

