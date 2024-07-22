#include "kalman.h"

namespace intnavlib {

Eigen::Matrix<double,15,15> InitializeLcPMmatrix(const LcKfConfig & lc_kf_config) {

    Eigen::Matrix<double,15,15> P_matrix;

    // Initialize error covariance matrix
    P_matrix = Eigen::Matrix<double,15,15>::Zero();
    P_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_att_unc,2); // attitude error
    P_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_vel_unc,2); // vel error
    P_matrix.block<3,3>(6,6) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_pos_unc,2); // pos error
    P_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_b_a_unc,2); // acc bias error
    P_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_b_g_unc,2); // gyro bias error

    return P_matrix;
}

// Propagates error-state uncertainty, given old P and integration time
Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const LcKfConfig & lc_kf_config,
                                        const double & tor_i) {

    double geocentric_radius = R_0 / sqrt(1.0 - pow(e * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - e*e, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Eigen::Vector3d omega_ie_vec;
    omega_ie_vec << 0.0,0.0,omega_ie;
    Eigen::Matrix3d Omega_ie = skewSymmetric(omega_ie_vec);

    // 1. Determine error-state transition matrix using (14.50) (first-order approx)
    Eigen::Matrix<double,15,15> Phi_matrix = Eigen::Matrix<double,15,15>::Identity();
    Phi_matrix.block<3,3>(0,0) -= Omega_ie * tor_i;
    Phi_matrix.block<3,3>(0,12) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.block<3,3>(3,0) = - tor_i * skewSymmetric(old_nav_est_ecef.C_b_e * imu_meas.f);
    Phi_matrix.block<3,3>(3,3) -= 2.0 * Omega_ie * tor_i;
    Phi_matrix.block<3,3>(3,6) = -tor_i * 2 * gravityEcef(old_nav_est_ecef.r_eb_e) /
        geocentric_radius * old_nav_est_ecef.r_eb_e.transpose() / sqrt(old_nav_est_ecef.r_eb_e.squaredNorm());
    Phi_matrix.block<3,3>(3,9) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.block<3,3>(6,3) = Eigen::Matrix3d::Identity() * tor_i;

    // 2. Determine approximate system noise covariance matrix using (14.82)
    Eigen::Matrix<double,15,15> Q_prime_matrix = Eigen::Matrix<double,15,15>::Zero();
    Q_prime_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * lc_kf_config.gyro_noise_PSD * tor_i;
    Q_prime_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * lc_kf_config.accel_noise_PSD * tor_i;
    Q_prime_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * lc_kf_config.accel_bias_PSD * tor_i;
    Q_prime_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * lc_kf_config.gyro_bias_PSD * tor_i;

    // 4. Propagate state estimation error covariance matrix using (3.46)
    Eigen::Matrix<double,15,15> P_matrix = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
    return P_matrix;
}

StateEstEcefLc lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                    const Eigen::Matrix<double,15,15> & P_matrix_propagated,
                                    const StateEstEcefLc & state_est_old) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_propagated = Eigen::Matrix<double,15,1>::Zero();
        
    // 5. Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,3,15> H_matrix = Eigen::Matrix<double,3,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position

    // 6. Set-up measurement noise covariance matrix
    Eigen::Matrix3d R_matrix = pos_meas.cov_mat;

    // 7. Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, 15,3> K_matrix = P_matrix_propagated * H_matrix.transpose() * 
                                (H_matrix * P_matrix_propagated * H_matrix.transpose() + R_matrix).inverse();

    // 8. Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Eigen::Vector3d delta_z = pos_meas.r_eb_e - state_est_old.nav_sol.r_eb_e;

    // 9. Update error state estimates using (3.24)
    Eigen::Matrix<double,15,1> x_est_new = x_est_propagated + K_matrix * delta_z;

    // 10. Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_new = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * P_matrix_propagated;

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcefLc state_est_new;
    state_est_new.nav_sol.time = state_est_old.nav_sol.time;
    state_est_new.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_old.nav_sol.C_b_e;
    state_est_new.nav_sol.v_eb_e = state_est_old.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_new.nav_sol.r_eb_e = state_est_old.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);

    // Update IMU bias estimates
    state_est_new.acc_bias += x_est_new.block<3,1>(9,0);
    state_est_new.gyro_bias += x_est_new.block<3,1>(12,0);

    return state_est_new;

}

};
