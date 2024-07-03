#include "helpers.h"

namespace helpers {

NavSolutionEcef nedToEcef(const NavSolutionNed& nav_sol_ned) {

    // Extract the inputs
    double L_b = nav_sol_ned.latitude;
    double lambda_b = nav_sol_ned.longitude;
    double h_b = nav_sol_ned.height;
    Eigen::Vector3d v_eb_n = nav_sol_ned.v_eb_n;
    Eigen::Matrix3d C_b_n = nav_sol_ned.C_b_n;

    // Compute transverse radius of curvature
    double sin_lat = std::sin(L_b);
    double cos_lat = std::cos(L_b);
    double sin_long = std::sin(lambda_b);
    double cos_long = std::cos(lambda_b);
    double R_E = R_0 / std::sqrt(1 - (e * sin_lat) * (e * sin_lat));

    // Convert position from curvilinear to Cartesian ECEF coordinates
    Eigen::Vector3d r_eb_e;
    r_eb_e(0) = (R_E + h_b) * cos_lat * cos_long;
    r_eb_e(1) = (R_E + h_b) * cos_lat * sin_long;
    r_eb_e(2) = ((1 - e * e) * R_E + h_b) * sin_lat;

    // Compute the ECEF to NED coordinate transformation matrix
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
             -sin_long,            cos_long,          0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;

    // Transform velocity from NED to ECEF frame
    Eigen::Vector3d v_eb_e = C_e_n.transpose() * v_eb_n;

    // Transform attitude from NED to ECEF frame
    Eigen::Matrix3d C_b_e = C_e_n.transpose() * C_b_n;

    // Construct the output structure
    NavSolutionEcef nav_sol_ecef;
    nav_sol_ecef.time = nav_sol_ned.time;
    nav_sol_ecef.p_eb_e = r_eb_e;
    nav_sol_ecef.v_eb_e = v_eb_e;
    nav_sol_ecef.C_b_e = C_b_e;

    return nav_sol_ecef;
}

NavSolutionNed ecefToNed(const NavSolutionEcef & nav_sol_ecef){

// Convert position using Borkowski closed-form exact solution
// From (2.113)
double lambda_b = atan2(nav_sol_ecef.p_eb_e(1), nav_sol_ecef.p_eb_e(0));

// From (C.29) and (C.30)
double k1 = sqrt(1 - pow(e,2)) * abs(nav_sol_ecef.p_eb_e(2));
double k2 = pow(e,2) * R_0;
double beta = sqrt(pow(nav_sol_ecef.p_eb_e(0),2) + pow(nav_sol_ecef.p_eb_e(1),2));
double E = (k1 - k2) / beta;
double F = (k1 + k2) / beta;

// From (C.31)
double P = 4/3 * (E*F + 1);

// From (C.32)
double Q = 2 * (pow(E,2) - pow(F,2));

// From (C.33)
double D = pow(P,3) + pow(Q,2);

// From (C.34)
double V = pow(sqrt(D) - Q, 1/3) - pow(sqrt(D) + Q,1/3);

// From (C.35)
double G = 0.5 * (sqrt(pow(E,2) + V) + E);

// From (C.36)
double T = sqrt(pow(G,2) + (F - V * G) / (2 * G - E)) - G;

// From (C.37)
double L_b = sgn(nav_sol_ecef.p_eb_e(2)) * atan((1 - pow(T,2)) / (2 * T * sqrt (1 - pow(e,2))));

// From (C.38)
double h_b = (beta - R_0 * T) * cos(L_b) +
    (nav_sol_ecef.p_eb_e(2) - sgn(nav_sol_ecef.p_eb_e(2)) * R_0 * sqrt(1 - pow(e,2))) * sin(L_b);
      
// Calculate ECEF to NED coordinate transformation matrix using (2.150)
double cos_lat = cos(L_b);
double sin_lat = sin(L_b);
double cos_long = cos(lambda_b);
double sin_long = sin(lambda_b);

Eigen::Matrix3d C_e_n;
C_e_n << -sin_lat * cos_long, -sin_lat * sin_long,  cos_lat,
                   -sin_long,            cos_long,        0,
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


Eigen::Matrix3d rpyToR(const Eigen::Vector3d & rpy) {

    double roll = rpy(0);
    double pitch = rpy(1);
    double yaw = rpy(2);

    // Precompute sines and cosines of Euler angles
    double sin_phi = std::sin(roll);
    double cos_phi = std::cos(roll);
    double sin_theta = std::sin(pitch);
    double cos_theta = std::cos(pitch);
    double sin_psi = std::sin(yaw);
    double cos_psi = std::cos(yaw);

    // Calculate the coordinate transformation matrix using the provided equations
    Eigen::Matrix3d C;
    C(0,0) = cos_theta * cos_psi;
    C(0,1) = cos_theta * sin_psi;
    C(0,2) = -sin_theta;
    C(1,0) = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
    C(1,1) = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
    C(1,2) = sin_phi * cos_theta;
    C(2,0) = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
    C(2,1) = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
    C(2,2) = cos_phi * cos_theta;

    return C;
}

Eigen::Vector3d rToRpy(const Eigen::Matrix3d & C) {

    Eigen::Vector3d rpy;

    rpy(0) = atan2(C(1,2),C(2,2));
    rpy(1) = - asin(C(0,2));      
    rpy(2) = atan2(C(0,1),C(0,0));

    return rpy;
}

ImuMeasurements kinematicsEcef(const NavSolutionEcef & old_nav, const NavSolutionEcef & new_nav) {
    
    double tor_i = new_nav.time - old_nav.time;

    // Init measurements to 0
    ImuMeasurements true_imu_meas;
    true_imu_meas.time = new_nav.time;
    true_imu_meas.f = Eigen::Vector3d::Zero();
    true_imu_meas.omega = Eigen::Vector3d::Zero();

    if (tor_i > 0) {

    // From (2.145) determine the Earth rotation over the update interval
    // C_Earth = C_e_i' * old_C_e_i

    double alpha_ie = omega_ie * tor_i;

    double cos_alpha_ie = cos(alpha_ie);
    double sin_alpha_ie = sin(alpha_ie);

    Eigen::Matrix3d C_Earth;
    C_Earth << cos_alpha_ie, sin_alpha_ie, 0,
                            -sin_alpha_ie, cos_alpha_ie, 0,
                                        0,             0,  1;

    // Obtain coordinate transformation matrix from the old attitude (w.r.t.
    // an inertial frame) to the new (compensate for earth rotation)
    Eigen::Matrix3d C_old_new = new_nav.C_b_e.transpose() * C_Earth * old_nav.C_b_e;

    // Calculate the approximate angular rate w.r.t. an inertial frame
    Eigen::Vector3d alpha_ib_b;
    alpha_ib_b(0) = 0.5 * (C_old_new(1,2) - C_old_new(2,1));
    alpha_ib_b(1) = 0.5 * (C_old_new(2,0) - C_old_new(0,2));
    alpha_ib_b(2) = 0.5 * (C_old_new(0,1) - C_old_new(1,0));

    // Calculate and apply the scaling factor
    double temp = acos(0.5 * (C_old_new(1,1) + C_old_new(2,2) + C_old_new(3,3) - 1.0));
    if (temp > 2e-5) // scaling is 1 if temp is less than this
        alpha_ib_b = alpha_ib_b * temp/sin(temp);
    
    // Calculate the angular rate
    true_imu_meas.omega = alpha_ib_b / tor_i;

    // Calculate the specific force resolved about ECEF-frame axes
    // From (5.36)
    Eigen::Vector3d omega_ie_vec;
    omega_ie_vec << 0 , 0 , omega_ie;
    Eigen::Vector3d f_ib_e = ((new_nav.v_eb_e - old_nav.v_eb_e) / tor_i) - gravityEcef(new_nav.p_eb_e)
        + 2 * skewSymmetric(omega_ie_vec) * old_nav.v_eb_e;

    // Calculate the average body-to-ECEF-frame coordinate transformation
    // matrix over the update interval using (5.84) and (5.85)

    double mag_alpha = alpha_ib_b.norm();
    Eigen::Matrix3d Alpha_ib_b = skewSymmetric(alpha_ib_b);

    Eigen::Vector3d alpha_ie_vec;
    alpha_ie_vec << 0 , 0 , alpha_ie;   

    Eigen::Matrix3d ave_C_b_e;
    if (mag_alpha>1e-8) {
        ave_C_b_e = old_nav.C_b_e * 
            (Eigen::Matrix3d::Identity() + 

            (1 - cos(mag_alpha) / pow(mag_alpha,2)) *
            Alpha_ib_b + 
            
            (1 - sin(mag_alpha) / mag_alpha) / pow(mag_alpha,2) * 
            Alpha_ib_b * Alpha_ib_b) - 

            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    else { // Approximate if angle small enough (sum not multiply)
        ave_C_b_e = old_nav.C_b_e -
            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    
    // Transform specific force to body-frame resolving axes using (5.81)
    // Groves inverts, but cant we just transpose?
    true_imu_meas.f = ave_C_b_e.transpose() * f_ib_e;

    }

    return true_imu_meas;

}

Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e) {

    double mag_r = r_eb_e.norm();

    Eigen::Vector3d g = Eigen::Vector3d::Zero();
    // If the input position is 0,0,0, produce a dummy output

    if (mag_r!=0)
    // Calculate gravitational acceleration using (2.142)
    {
        double z_scale = 5 * pow(r_eb_e(2) / mag_r,2);

        Eigen::Vector3d gamma_1;
        gamma_1 << (1 - z_scale) * r_eb_e(0), 
                    (1 - z_scale) * r_eb_e(1),
                    (3 - z_scale) * r_eb_e(2);

        Eigen::Vector3d gamma;
        gamma = -mu / pow(mag_r,3) *
                (r_eb_e + 1.5 * J_2 * 
                pow(R_0 / mag_r,2) * gamma_1);

        // Add centripetal acceleration using (2.133)
        g(0) = gamma(0) + pow(omega_ie,2) * r_eb_e(0);
        g(1) = gamma(1) + pow(omega_ie,2) * r_eb_e(1);
        g(2) = gamma(2);
    }

    return g;

}

Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d & a) {
    
    Eigen::Matrix3d S;
    
    S << 0, -a(2),  a(1),
      a(2),     0, -a(0),
     -a(1),  a(0),     0;

    return S;

}

ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                        const ImuErrors & imu_errors,
                        const double & tor_i,
                        std::mt19937 & gen) {
    
    ImuMeasurements imu_measurements;

    imu_measurements.time = true_imu_meas.time;

    // Init noises
    Eigen::Vector3d accel_noise = Eigen::Vector3d::Zero();
    Eigen::Vector3d gyro_noise = Eigen::Vector3d::Zero();

    // Normal distribution init
    std::normal_distribution<double> randn(0, 1); 

    if(tor_i > 0) {
        accel_noise = Eigen::Vector3d::Ones() * randn(gen) * imu_errors.accel_noise_root_PSD / sqrt(tor_i);  
        gyro_noise = Eigen::Vector3d::Ones() * randn(gen) * imu_errors.gyro_noise_root_PSD / sqrt(tor_i);
    }

    // Calculate accelerometer and gyro outputs using (4.16) and (4.17)

    // Specific force
    Eigen::Vector3d uq_f_ib_b = imu_errors.b_a + (Eigen::Matrix3d::Identity() + imu_errors.M_a) * true_imu_meas.f 
                                + accel_noise;

    // Angular velocity
    Eigen::Vector3d uq_omega_ib_b = imu_errors.b_g + (Eigen::Matrix3d::Identity() + imu_errors.M_g) * true_imu_meas.omega 
                                    + imu_errors.G_g * true_imu_meas.f 
                                    + gyro_noise;

    // // Quantize accelerometer outputs
    // if (IMU_errors.accel_quant_level > 0) {
    //     meas_f_ib_b = IMU_errors.accel_quant_level * round((uq_f_ib_b +...
    //         old_quant_residuals_f) / IMU_errors.accel_quant_level);
    //     quant_residuals_f = uq_f_ib_b + old_quant_residuals_f -...
    //         meas_f_ib_b;
    // }
    // else{
    //     imu_measurements.f = uq_f_ib_b;
    //     quant_residuals_f = [0;0;0];
    // }

    // // Quantize gyro outputs
    // if (IMU_errors.gyro_quant_level>0) {
    //     imu_measurements.omega = IMU_errors.gyro_quant_level * round((uq_omega_ib_b +...
    //         old_quant_residuals_omega) / IMU_errors.gyro_quant_level);
    //     quant_residuals_omega = uq_omega_ib_b + old_quant_residuals_omega -...
    //         imu_measurements.omega;
    // }
    // else {
    //     imu_measurements.omega = uq_omega_ib_b;
    //     quant_residuals_omega = [0;0;0];
    // }

    // No quantization
    imu_measurements.f = uq_f_ib_b;
    imu_measurements.omega = uq_omega_ib_b;

    return imu_measurements;
}

NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i) {

NavSolutionEcef new_nav;

new_nav.time = imu_meas.time;

// ATTITUDE UPDATE
// From (2.145) determine the Earth rotation over the update interval
// C_Earth = C_e_i' * old_C_e_i

double alpha_ie = omega_ie * tor_i;

double cos_alpha_ie = cos(alpha_ie);
double sin_alpha_ie = sin(alpha_ie);

Eigen::Matrix3d C_Earth;
C_Earth << cos_alpha_ie, sin_alpha_ie, 0,
            -sin_alpha_ie, cos_alpha_ie, 0,
                        0,             0,  1;
                       
// Calculate attitude increment, magnitude, and skew-symmetric matrix
Eigen::Vector3d alpha_ib_b = imu_meas.omega * tor_i;
double mag_alpha = alpha_ib_b.norm();
Eigen::Matrix3d Alpha_ib_b = skewSymmetric(alpha_ib_b);  

// Obtain coordinate transformation matrix from the new attitude w.r.t. an
// inertial frame to the old using Rodrigues' formula, (5.73)
Eigen::Matrix3d C_new_old;
if (mag_alpha>1e-8) {
    C_new_old = Eigen::Matrix3d::Identity() + 
                sin(mag_alpha) / mag_alpha * Alpha_ib_b +
                (1 - cos(mag_alpha)) / pow(mag_alpha,2) * Alpha_ib_b * Alpha_ib_b;
}
else {
    C_new_old = Eigen::Matrix3d::Identity() + Alpha_ib_b;
}

// Update attitude using (5.75)
new_nav.C_b_e = C_Earth * old_nav.C_b_e * C_new_old;
    
// SPECIFIC FORCE FRAME TRANSFORMATION
// Calculate the average body-to-ECEF-frame coordinate transformation
// matrix over the update interval using (5.84) and (5.85)

Eigen::Vector3d alpha_ie_vec;
alpha_ie_vec << 0 , 0 , alpha_ie;   

Eigen::Matrix3d ave_C_b_e;
if (mag_alpha>1e-8) {
    ave_C_b_e = old_nav.C_b_e * 
        (Eigen::Matrix3d::Identity() + 

        (1 - cos(mag_alpha) / pow(mag_alpha,2)) *
        Alpha_ib_b + 
        
        (1 - sin(mag_alpha) / mag_alpha) / pow(mag_alpha,2) * 
        Alpha_ib_b * Alpha_ib_b) - 

        0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
}
else { // Approximate if angle small enough (sum not multiply)
    ave_C_b_e = old_nav.C_b_e -
        0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
}

// Transform specific force to ECEF-frame resolving axes using (5.85)
Eigen::Vector3d f_ib_e = ave_C_b_e * imu_meas.f;
    
// UPDATE VELOCITY
// From (5.36)

Eigen::Vector3d omega_ie_vec;
alpha_ie_vec << 0 , 0 , omega_ie;   

new_nav.v_eb_e = old_nav.v_eb_e + tor_i * (f_ib_e + gravityEcef(old_nav.p_eb_e) -
                        2 * skewSymmetric(omega_ie_vec) * old_nav.v_eb_e);

// UPDATE CARTESIAN POSITION
// From (5.38),
new_nav.p_eb_e = old_nav.p_eb_e + (new_nav.v_eb_e + old_nav.v_eb_e) * 0.5 * tor_i; 

return new_nav;

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
    double temp = 1 - pow((e * sin(L)),2); 
    double R_N = R_0 * (1 - pow(e,2)) / pow(temp,1.5);
    // Calculate transverse radius of curvature using (2.105)
    double R_E = R_0 / sqrt(temp);
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
    Eigen::Vector3d delta_eul_nb_n = -rToRpy(delta_C_b_n);

    ErrorsNed errors_ned;
    errors_ned.time = true_nav_sol.time;
    errors_ned.delta_r_eb_n = delta_r_eb_n;
    errors_ned.delta_v_eb_n = delta_v_eb_n;
    errors_ned.delta_eul_nb_n = delta_eul_nb_n;

    return errors_ned;
}

};