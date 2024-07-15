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
                                        const double & tor_s) {

    double geocentric_radius = R_0 / sqrt(1.0 - pow(e * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - e*e, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Eigen::Vector3d omega_ie_vec;
    omega_ie_vec << 0.0,0.0,omega_ie;
    Eigen::Matrix3d Omega_ie = skewSymmetric(omega_ie_vec);

    // 1. Determine error-state transition matrix using (14.50) (first-order approx)
    Eigen::Matrix<double,15,15> Phi_matrix = Eigen::Matrix<double,15,15>::Identity();
    Phi_matrix.block<3,3>(0,0) -= Omega_ie * tor_s;
    Phi_matrix.block<3,3>(0,12) = old_nav_est_ecef.C_b_e * tor_s;
    Phi_matrix.block<3,3>(3,12) -= tor_s * skewSymmetric(old_nav_est_ecef.C_b_e * imu_meas.f);
    Phi_matrix.block<3,3>(3,3) -= 2.0 * Omega_ie * tor_s;
    Phi_matrix.block<3,3>(3,6) = -tor_s * 2 * gravityEcef(old_nav_est_ecef.r_eb_e) /
        geocentric_radius * old_nav_est_ecef.r_eb_e.transpose() / sqrt(old_nav_est_ecef.r_eb_e.norm());
    Phi_matrix.block<3,3>(3,9) = old_nav_est_ecef.C_b_e * tor_s;
    Phi_matrix.block<3,3>(6,3) = Eigen::Matrix3d::Identity() * tor_s;

    // 2. Determine approximate system noise covariance matrix using (14.82)
    Eigen::Matrix<double,15,15> Q_prime_matrix = Eigen::Matrix<double,15,15>::Identity();
    Q_prime_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * lc_kf_config.gyro_noise_PSD * tor_s;
    Q_prime_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * lc_kf_config.accel_noise_PSD * tor_s;
    Q_prime_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * lc_kf_config.accel_bias_PSD * tor_s;
    Q_prime_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * lc_kf_config.gyro_bias_PSD * tor_s;

    // 4. Propagate state estimation error covariance matrix using (3.46)
    Eigen::Matrix<double,15,15> P_matrix = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
    return P_matrix;
}

};
