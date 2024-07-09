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

    double alpha_ie = omega_ie * tor_i;

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
    Eigen::Vector3d omega_ie_vec;
    omega_ie_vec << 0.0 , 0.0 , omega_ie;
    Eigen::Vector3d f_ib_e = ((new_nav.v_eb_e - old_nav.v_eb_e) / tor_i) - gravityEcef(old_nav.r_eb_e)
        + 2.0 * skewSymmetric(omega_ie_vec) * old_nav.v_eb_e;

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
    true_imu_meas.f = ave_C_b_e.inverse() * f_ib_e; // transpose or inverse?

    }

    return true_imu_meas;

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
    std::normal_distribution<double> randn(0.0, 1.0); 

    if(tor_i > 0.0) {
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


};

