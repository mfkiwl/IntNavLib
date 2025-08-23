#include "helpers.h"

namespace intnavlib {

NavSolutionEcef nedToEcef(const NavSolutionNed& nav_sol_ned) {
    // Extract the inputs
    double L_b = nav_sol_ned.latitude;
    double lambda_b = nav_sol_ned.longitude;
    double h_b = nav_sol_ned.height;
    Eigen::Vector3d v_eb_n = nav_sol_ned.v_eb_n;
    Eigen::Matrix3d C_b_n = nav_sol_ned.C_b_n;
    double sin_lat = std::sin(L_b);
    double cos_lat = std::cos(L_b);
    double sin_long = std::sin(lambda_b);
    double cos_long = std::cos(lambda_b);
    // Compute transverse radius of curvature
    double R_E = kR0 / std::sqrt(1.0 - pow((kEccentricity * sin_lat),2.0));
    // Convert position from curvilinear to Cartesian ECEF coordinates
    Eigen::Vector3d r_eb_e;
    r_eb_e(0) = (R_E + h_b) * cos_lat * cos_long;
    r_eb_e(1) = (R_E + h_b) * cos_lat * sin_long;
    r_eb_e(2) = ((1.0 - kEccentricity * kEccentricity) * R_E + h_b) * sin_lat;
    // Compute the ECEF to NED coordinate transformation matrix
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
             -sin_long,            cos_long,          0.0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    // Transform velocity from NED to ECEF frame
    Eigen::Vector3d v_eb_e = C_e_n.transpose() * v_eb_n;
    // Transform attitude from NED to ECEF frame
    Eigen::Matrix3d C_b_e = C_e_n.transpose() * C_b_n;
    // Construct the output structure
    NavSolutionEcef nav_sol_ecef;
    nav_sol_ecef.time = nav_sol_ned.time;
    nav_sol_ecef.r_eb_e = r_eb_e;
    nav_sol_ecef.v_eb_e = v_eb_e;
    nav_sol_ecef.C_b_e = C_b_e;
    return nav_sol_ecef;
}

NavSolutionNed ecefToNed(const NavSolutionEcef & nav_sol_ecef){
    // Convert position using Borkowski closed-form exact solution
    // From (2.113)
    double lambda_b = atan2(nav_sol_ecef.r_eb_e(1), nav_sol_ecef.r_eb_e(0));
    // From (C.29) and (C.30)
    double k1 = sqrt(1.0 - pow(kEccentricity,2.0)) * abs(nav_sol_ecef.r_eb_e(2));
    double k2 = pow(kEccentricity,2.0) * kR0;
    double beta = sqrt(pow(nav_sol_ecef.r_eb_e(0),2.0) + pow(nav_sol_ecef.r_eb_e(1),2.0));
    double E = (k1 - k2) / beta;
    double F = (k1 + k2) / beta;
    // From (C.31)
    double P = 4.0/3.0 * (E*F + 1.0);
    // From (C.32)
    double Q = 2.0 * (pow(E,2.0) - pow(F,2.0));
    // From (C.33)
    double D = pow(P,3.0) + pow(Q,2.0);
    // From (C.34)
    double V = pow(sqrt(D) - Q, 1.0/3.0) - pow(sqrt(D) + Q,1.0/3.0);
    // From (C.35)
    double G = 0.5 * (sqrt(pow(E,2.0) + V) + E);
    // From (C.36)
    double T = sqrt(pow(G,2.0) + (F - V * G) / (2.0 * G - E)) - G;
    // From (C.37)
    double L_b = sgn(nav_sol_ecef.r_eb_e(2)) * atan((1.0 - pow(T,2.0)) / (2.0 * T * sqrt (1.0 - pow(kEccentricity,2.0))));
    // From (C.38)
    double h_b = (beta - kR0 * T) * cos(L_b) +
        (nav_sol_ecef.r_eb_e(2) - sgn(nav_sol_ecef.r_eb_e(2)) * kR0 * sqrt(1.0 - pow(kEccentricity,2.0))) * sin(L_b);  
    // Calculate ECEF to NED coordinate transformation matrix using (2.150)
    double cos_lat = cos(L_b);
    double sin_lat = sin(L_b);
    double cos_long = cos(lambda_b);
    double sin_long = sin(lambda_b);
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
                    -sin_long,            cos_long,        0.0,
            -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;
    // Transform velocity using (2.73)
    Eigen::Vector3d v_eb_n = C_e_n * nav_sol_ecef.v_eb_e;
    // Transform attitude using (2.15)
    Eigen::Matrix3d C_b_n = C_e_n * nav_sol_ecef.C_b_e;
    NavSolutionNed nav_sol_ned;
    nav_sol_ned.time = nav_sol_ecef.time;
    nav_sol_ned.latitude = L_b;
    nav_sol_ned.longitude = lambda_b;
    nav_sol_ned.height = h_b;
    nav_sol_ned.v_eb_n = v_eb_n;
    nav_sol_ned.C_b_n = C_b_n;
    return nav_sol_ned;
}

ErrorsSigmasEcef getErrorsSigmasEcef(const StateEstEcef & state_est_ecef, const NavSolutionEcef & true_nav_ecef) {   
    ErrorsSigmasEcef errors_sigmas;
    errors_sigmas.time = state_est_ecef.nav_sol.time;
    errors_sigmas.delta_r_eb_e = state_est_ecef.nav_sol.r_eb_e - true_nav_ecef.r_eb_e;
    errors_sigmas.delta_v_eb_e = state_est_ecef.nav_sol.v_eb_e - true_nav_ecef.v_eb_e;
    errors_sigmas.delta_rot_eb_e = deSkew(state_est_ecef.nav_sol.C_b_e * true_nav_ecef.C_b_e.transpose() - Eigen::Matrix3d::Identity());
    
    errors_sigmas.sigma_delta_r_eb_e << sqrt(state_est_ecef.P_matrix(6,6)), sqrt(state_est_ecef.P_matrix(7,7)), sqrt(state_est_ecef.P_matrix(8,8));
    errors_sigmas.sigma_delta_v_eb_e << sqrt(state_est_ecef.P_matrix(3,3)), sqrt(state_est_ecef.P_matrix(4,4)), sqrt(state_est_ecef.P_matrix(5,5));
    errors_sigmas.sigma_delta_rot_eb_e << sqrt(state_est_ecef.P_matrix(0,0)), sqrt(state_est_ecef.P_matrix(1,1)), sqrt(state_est_ecef.P_matrix(2,2));
    return errors_sigmas;
}

Eigen::Matrix3d eulerToDcm(const Eigen::Vector3d & rpy) {
    double roll = rpy(0);
    double pitch = rpy(1);
    double yaw = rpy(2);
    // Precompute sines and cosines of Euler angles
    double sin_roll = std::sin(roll);
    double cos_roll = std::cos(roll);
    double sin_pitch = std::sin(pitch);
    double cos_pitch = std::cos(pitch);
    double sin_yaw = std::sin(yaw);
    double cos_yaw = std::cos(yaw);
    // Calculate the coordinate transformation matrix R = Rz(yaw) * Ry(pitch) * Rx(roll)
    Eigen::Matrix3d C;
    C(0,0) = cos_pitch * cos_yaw;
    C(1,0) = cos_pitch * sin_yaw;
    C(2,0) = -sin_pitch;
    C(0,1) = -cos_roll * sin_yaw + sin_roll * sin_pitch * cos_yaw;
    C(1,1) = cos_roll * cos_yaw + sin_roll * sin_pitch * sin_yaw;
    C(2,1) = sin_roll * cos_pitch;
    C(0,2) = sin_roll * sin_yaw + cos_roll * sin_pitch * cos_yaw;
    C(1,2) = -sin_roll * cos_yaw + cos_roll * sin_pitch * sin_yaw;
    C(2,2) = cos_roll * cos_pitch;
    return C;
}

Eigen::Vector3d dcmToEuler(const Eigen::Matrix3d & C) {
    Eigen::Vector3d rpy;
    rpy(0) = atan2(C(2,1),C(2,2));
    rpy(1) = - asin(C(2,0));      
    rpy(2) = atan2(C(1,0),C(0,0));
    return rpy;
}

Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e) {
    double mag_r = r_eb_e.norm();
    Eigen::Vector3d g = Eigen::Vector3d::Zero();
    // If the input position is 0,0,0, produce a dummy output
    if (mag_r >= kEpsilon)
    // Calculate gravitational acceleration using (2.142)
    {
        double z_scale = 5.0 * pow(r_eb_e(2) / mag_r,2.0);
        Eigen::Vector3d gamma_1;
        gamma_1 << (1.0 - z_scale) * r_eb_e(0), 
                    (1.0 - z_scale) * r_eb_e(1),
                    (3.0 - z_scale) * r_eb_e(2);
        Eigen::Vector3d gamma;
        gamma = (-kGravConst / pow(mag_r,3.0)) *
                (r_eb_e + 1.5 * kJ2 * 
                pow(kR0 / mag_r,2.0) * gamma_1);
        // Add centripetal acceleration using (2.133)
        g(0) = gamma(0) + pow(kOmega_ie,2.0) * r_eb_e(0);
        g(1) = gamma(1) + pow(kOmega_ie,2.0) * r_eb_e(1);
        g(2) = gamma(2);
    }
    return g;
}

Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d & a) {
    Eigen::Matrix3d S;
    S << 0.0, -a(2),  a(1),
      a(2),     0.0, -a(0),
     -a(1),  a(0),     0.0;
    return S;
}

Eigen::Vector3d deSkew(const Eigen::Matrix3d & S){
    Eigen::Vector3d a;
    a << -S(1,2), S(0,2), -S(0,1);
    return a; 
}

// Function to get the current date and time as a formatted string
std::string getCurrentDateTime() {
    auto now = std::time(nullptr);
    std::tm tm_now;
    localtime_r(&now, &tm_now); // Use localtime_s on Windows or localtime_r on Unix-like systems
    std::ostringstream oss;
    oss << std::put_time(&tm_now, "%Y%m%d_%H%M%S");
    return oss.str();
}

Eigen::Vector2d radiiOfCurvature(double L) {
    // Calculate meridian radius of curvature using (2.105)
    double temp = 1.0 - pow((kEccentricity * sin(L)),2.0); 
    double R_N = kR0 * (1.0 - pow(kEccentricity,2.0)) / pow(temp,1.5);
    // Calculate transverse radius of curvature using (2.105)
    double R_E = kR0 / sqrt(temp);
    Eigen::Vector2d radii;
    radii << R_N, R_E;
    return radii;
}

ErrorsNed calculateErrorsNed(const NavSolutionNed & true_nav_sol, 
                                const NavSolutionNed & est_nav_sol){
    // Position error calculation
    Eigen::Vector2d radii = radiiOfCurvature(true_nav_sol.latitude);
    double R_N = radii(0);
    double R_E = radii(1);
    Eigen::Vector3d delta_r_eb_n;
    delta_r_eb_n(0) = (est_nav_sol.latitude - true_nav_sol.latitude) * (R_N + true_nav_sol.height);
    delta_r_eb_n(1) = (est_nav_sol.longitude - true_nav_sol.longitude) * (R_E + true_nav_sol.height) * std::cos(true_nav_sol.latitude);
    delta_r_eb_n(2) = -(est_nav_sol.height - true_nav_sol.height);
    // Velocity error calculation
    Eigen::Vector3d delta_v_eb_n = est_nav_sol.v_eb_n - true_nav_sol.v_eb_n;
    // Attitude error calculation
    Eigen::Matrix3d delta_C_b_n = est_nav_sol.C_b_n * true_nav_sol.C_b_n.transpose();
    Eigen::Vector3d delta_rot_nb_n = -dcmToEuler(delta_C_b_n);
    ErrorsNed errors_ned;
    errors_ned.time = true_nav_sol.time;
    errors_ned.delta_r_eb_n = delta_r_eb_n;
    errors_ned.delta_v_eb_n = delta_v_eb_n;
    errors_ned.delta_rot_nb_n = delta_rot_nb_n;
    return errors_ned;
}

ImuErrors tacticalImuErrors(){
    ImuErrors imu_errors;
    imu_errors.b_a << 900.0,-1300.0,800.0;
    imu_errors.b_a = imu_errors.b_a * kMuGToMetersPerSecondSquared;
    imu_errors.b_g << -9.0, 13.0, -8.0;
    imu_errors.b_g = imu_errors.b_g * kDegToRad / 3600.0;
    imu_errors.M_a << 500.0, -300.0, 200.0,
            -150.0, -600.0, 250.0,
            -250.0,  100.0, 450.0;
    imu_errors.M_a = imu_errors.M_a * 1.0e-6;
    imu_errors.M_g << 400.0, -300.0,  250.0,
            0.0, -300.0, -150.0,
            0.0,    0.0, -350.0; 
    imu_errors.M_g = imu_errors.M_g * 1.0e-6;
    imu_errors.G_g << 0.9, -1.1, -0.6,
            -0.5,  1.9, -1.6,
            0.3,  1.1, -1.3;
    imu_errors.G_g = imu_errors.G_g * kDegToRad / (3600.0 * 9.80665);  
    imu_errors.accel_noise_root_psd = 100.0 * kMuGToMetersPerSecondSquared;
    imu_errors.gyro_noise_root_psd = 0.01 * kDegToRad / 60.0;
    imu_errors.accel_quant_level = 1.0e-2;
    imu_errors.gyro_quant_level = 2.0e-4;
    return imu_errors;
}

GnssConfig defaultGnssConfig(){
    GnssConfig gnss_config;
    gnss_config.epoch_interval = 0.5;
    gnss_config.init_est_r_ea_e = Eigen::Vector3d::Zero();
    gnss_config.no_sat = 30.0;
    gnss_config.r_os = 2.656175E7;
    gnss_config.inclination = 55.0;
    gnss_config.const_delta_lambda = 0.0;
    gnss_config.const_delta_t = 0.0;
    gnss_config.mask_angle = 10.0;
    gnss_config.sis_err_sd = 1.0;
    gnss_config.zenith_iono_err_sd = 2.0;
    gnss_config.zenith_trop_err_sd = 0.2;
    gnss_config.code_track_err_sd = 1.0;
    gnss_config.rate_track_err_sd = 0.02;
    gnss_config.rx_clock_offset = 10000.0;
    gnss_config.rx_clock_drift = 100.0;
    gnss_config.lc_pos_sd = 2.5;
    gnss_config.lc_vel_sd = 0.1;
    gnss_config.pseudo_range_sd = 2.5;
    gnss_config.range_rate_sd = 0.1;
    return gnss_config;
}

KfConfig tacticalImuKFConfig(){
    KfConfig kf_config;
    kf_config.init_att_unc = kDegToRad * 1.0;
    kf_config.init_vel_unc = 0.1;
    kf_config.init_pos_unc = 10.0;
    kf_config.init_b_a_unc = 1000.0 * kMuGToMetersPerSecondSquared;
    kf_config.init_b_g_unc = 10.0 * kDegToRad / 3600.0;
    kf_config.init_clock_offset_unc = 10.0;
    kf_config.init_clock_drift_unc = 0.1;
    kf_config.gyro_noise_psd = pow(0.02 * kDegToRad / 60.0, 2.0);
    kf_config.accel_noise_psd = pow(200.0 * kMuGToMetersPerSecondSquared, 2.0);
    kf_config.accel_bias_psd = 1.0E-7;
    kf_config.gyro_bias_psd = 2.0E-12;
    kf_config.clock_freq_psd = 1;
    kf_config.clock_phase_psd = 1;
    return kf_config;
}     


};