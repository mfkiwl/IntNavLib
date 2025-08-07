#include "kalman.h"

namespace intnavlib {

Eigen::Matrix<double,15,15> initializeLcPMmatrix(const KfConfig & lc_kf_config) {

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

Eigen::Matrix<double,17,17> initializeTcPMmatrix(const KfConfig & tc_kf_config) {

    Eigen::Matrix<double,17,17> P_matrix;

    // Initialize error covariance matrix
    P_matrix = Eigen::Matrix<double,17,17>::Zero();
    P_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * pow(tc_kf_config.init_att_unc,2); // attitude error
    P_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * pow(tc_kf_config.init_vel_unc,2); // vel error
    P_matrix.block<3,3>(6,6) = Eigen::Matrix3d::Identity() * pow(tc_kf_config.init_pos_unc,2); // pos error
    P_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * pow(tc_kf_config.init_b_a_unc,2); // acc bias error
    P_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * pow(tc_kf_config.init_b_g_unc,2); // gyro bias error
    P_matrix(15,15) = pow(tc_kf_config.init_clock_offset_unc,2); // clock offset error
    P_matrix(16,16) = pow(tc_kf_config.init_clock_drift_unc,2); // clock drift error

    return P_matrix;
}

// Propagates error-state uncertainty, given old P and integration time
Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & lc_kf_config,
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
        geocentric_radius * old_nav_est_ecef.r_eb_e.transpose() / old_nav_est_ecef.r_eb_e.norm();
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

Eigen::Matrix<double,17,17> tcPropUnc(const Eigen::Matrix<double,17,17> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & tc_kf_config,
                                        const double & tor_i) {
    double geocentric_radius = R_0 / sqrt(1.0 - pow(e * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - e*e, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Eigen::Vector3d omega_ie_vec;
    omega_ie_vec << 0.0,0.0,omega_ie;
    Eigen::Matrix3d Omega_ie = skewSymmetric(omega_ie_vec);

    // 1. Determine error-state transition matrix using (14.50) (first-order approx)
    Eigen::Matrix<double,17,17> Phi_matrix = Eigen::Matrix<double,17,17>::Identity();
    Phi_matrix.block<3,3>(0,0) -= Omega_ie * tor_i;
    Phi_matrix.block<3,3>(0,12) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.block<3,3>(3,0) = - tor_i * skewSymmetric(old_nav_est_ecef.C_b_e * imu_meas.f);
    Phi_matrix.block<3,3>(3,3) -= 2.0 * Omega_ie * tor_i;
    Phi_matrix.block<3,3>(3,6) = -tor_i * 2 * gravityEcef(old_nav_est_ecef.r_eb_e) /
        geocentric_radius * old_nav_est_ecef.r_eb_e.transpose() / old_nav_est_ecef.r_eb_e.norm();
    Phi_matrix.block<3,3>(3,9) = old_nav_est_ecef.C_b_e * tor_i;
    Phi_matrix.block<3,3>(6,3) = Eigen::Matrix3d::Identity() * tor_i;

    // Clock offset
    Phi_matrix(15,16) = tor_i;

    // 2. Determine approximate system noise covariance matrix using (14.82)
    Eigen::Matrix<double,17,17> Q_prime_matrix = Eigen::Matrix<double,17,17>::Zero();
    Q_prime_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * tc_kf_config.gyro_noise_PSD * tor_i;
    Q_prime_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * tc_kf_config.accel_noise_PSD * tor_i;
    Q_prime_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * tc_kf_config.accel_bias_PSD * tor_i;
    Q_prime_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * tc_kf_config.gyro_bias_PSD * tor_i;

    // Clock offset, drift
    Q_prime_matrix(15,15) = tc_kf_config.clock_phase_PSD * tor_i;
    Q_prime_matrix(16,16) = tc_kf_config.clock_freq_PSD * tor_i;

    // 4. Propagate state estimation error covariance matrix using (3.46)
    Eigen::Matrix<double,17,17> P_matrix = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
    return P_matrix;
}

StateEstEcefLc lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                const StateEstEcefLc & state_est_prior,
                                const double & p_value) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_propagated = Eigen::Matrix<double,15,1>::Zero();
        
    // Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,3,15> H_matrix = Eigen::Matrix<double,3,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position

    // Set-up measurement noise covariance matrix
    Eigen::Matrix3d R_matrix = pos_meas.cov_mat;

    // Calculate Kalman gain using (3.21)
    Eigen::Matrix3d S_matrix_inv = (H_matrix * state_est_prior.P_matrix * H_matrix.transpose() + R_matrix).inverse();
    Eigen::Matrix<double, 15,3> K_matrix = state_est_prior.P_matrix * H_matrix.transpose() * S_matrix_inv;

    // Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Eigen::Vector3d delta_z = pos_meas.r_eb_e - state_est_prior.nav_sol.r_eb_e;

    // Update error state estimates using (3.24)
    Eigen::Matrix<double,15,1> x_est_new = x_est_propagated + K_matrix * delta_z;

    // Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_post = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix;

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcefLc state_est_post;
    state_est_post.valid = true;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    // Real-time consistency check
    // See: Estimation with Applications to Tracking and Navigation -- Yaakov Bar-Shalom et al. p.237
    // The quadratic form res' * inv(S) * res is a determination of the standard Chi squared distribution with one dof.
    // We check it's inside the desired p_value
    // Simplified version: no history, just current measurement, 1 * n_res dof's.
    double chi2 = delta_z.dot(S_matrix_inv * delta_z);
    boost::math::chi_squared chi2_dist(3);
    double chi2_thresh = boost::math::quantile(chi2_dist, p_value);
    if (chi2 > chi2_thresh)
        // Set invalid update flag
        state_est_post.valid = false;
    
    // Return 
    return state_est_post;
}

StateEstEcefLc lcUpdateKFGnssEcef (const GnssPosVelMeasEcef & pos_vel_gnss_meas, 
                                    const StateEstEcefLc & state_est_prior) {
        
    // 5. Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,6,15> H_matrix = Eigen::Matrix<double,6,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position
    H_matrix.block<3,3>(3,3) = - Eigen::Matrix3d::Identity(); // Velocity

    // 6. Set-up measurement noise covariance matrix
    Eigen::Matrix<double,6,6> R_matrix = pos_vel_gnss_meas.cov_mat;

    // 7. Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, 15,6> K_matrix = state_est_prior.P_matrix * H_matrix.transpose() * 
                                (H_matrix * state_est_prior.P_matrix * H_matrix.transpose() + R_matrix).inverse();

    // 8. Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Eigen::Matrix<double,6,1> delta_z;
    delta_z.block<3,1>(0,0) = pos_vel_gnss_meas.r_ea_e - state_est_prior.nav_sol.r_eb_e;
    delta_z.block<3,1>(3,0) = pos_vel_gnss_meas.v_ea_e - state_est_prior.nav_sol.v_eb_e;

    // 9. Update error state estimates using (3.24)
    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_new = /*x_est_propagated + */ K_matrix * delta_z;

    // 10. Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_post = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix;

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcefLc state_est_post;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    return state_est_post;
}


StateEstEcefTc tcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcefTc & state_est_prior,
                                    const double & tor_s) {

    // Compute predicted range and range rate measurements from prior state est
    Eigen::Matrix<double, Eigen::Dynamic, 2, 0, MAX_GNSS_SATELLITES, 2> pred_meas(gnss_meas.no_meas,2);
    // Line of sight unit vectors in ecef
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, MAX_GNSS_SATELLITES, 3> u_as_e_T(gnss_meas.no_meas,3);

    // For each satellite
    for(int j=0; j<gnss_meas.no_meas; j++) {

        // Predicted range, approximated
        Eigen::Vector3d delta_r = gnss_meas.meas.block<1,3>(j,2).transpose() - state_est_prior.nav_sol.r_eb_e;
        double range = delta_r.norm();

        // Calculate ecef rotation during signal transit time using (8.36)
        Eigen::Matrix3d C_e_I;
        C_e_I << 1, omega_ie * range / c, 0,
                -omega_ie * range / c, 1, 0,
                0, 0, 1;

        // Predicted range, corrected
        // Convert satellite position in ecef frame at signal reception time, 
        // by taking into account earth rotation. This way you get the true range
        delta_r = C_e_I *  gnss_meas.meas.block<1,3>(j,2).transpose() - state_est_prior.nav_sol.r_eb_e;
        range = delta_r.norm();
        // Also consider prior on clock offset (meters) 
        pred_meas(j,0) = range + state_est_prior.clock_offset + state_est_prior.clock_drift * tor_s;

        // Predict pseudo-range rate using (9.165)
        // As before, get satellite velocity at signal reception time
        // Skew symmetric matrix of Earth rate
        Eigen::Vector3d omega_ie_vec;
        omega_ie_vec << 0.0,0.0,omega_ie;
        Eigen::Matrix3d Omega_ie = skewSymmetric(omega_ie_vec);
        // Get line of sight unit vector
        u_as_e_T.block<1,3>(j,0) = delta_r / range;
        double range_rate = u_as_e_T.block<1,3>(j,0) * (C_e_I * (gnss_meas.meas.block<1,3>(j,5).transpose() +
                                                        Omega_ie * gnss_meas.meas.block<1,3>(j,2).transpose()) - 

                                                        (state_est_prior.nav_sol.v_eb_e + 
                                                        Omega_ie * state_est_prior.nav_sol.r_eb_e)); 
        // Also consider prior on clock drift (meters/seconds) 
        pred_meas(j,1) = range_rate + state_est_prior.clock_drift;
    }

    // 5. Set-up measurement matrix using (14.126) - simplified Jacobian
    Eigen::Matrix<double, Eigen::Dynamic, 17, 0, 2 * MAX_GNSS_SATELLITES, 17> H_matrix = 
        Eigen::Matrix<double, Eigen::Dynamic, 17, 0, 2 * MAX_GNSS_SATELLITES, 17>::Zero(2*gnss_meas.no_meas, 17);

    // Ranges
    H_matrix.block(0,6,gnss_meas.no_meas,3) = u_as_e_T.block(0,0,gnss_meas.no_meas,3);
    H_matrix.block(0,15,gnss_meas.no_meas,1) = Eigen::MatrixXd::Ones(gnss_meas.no_meas,1);
    // Range rates
    H_matrix.block(gnss_meas.no_meas,3, gnss_meas.no_meas,3) = u_as_e_T.block(0,0,gnss_meas.no_meas, 3);
    H_matrix.block(gnss_meas.no_meas,16, gnss_meas.no_meas,1) = Eigen::MatrixXd::Ones(gnss_meas.no_meas,1);

    
    // 6. Set-up measurement noise covariance matrix
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* MAX_GNSS_SATELLITES, 2* MAX_GNSS_SATELLITES> R_matrix = gnss_meas.cov_mat;

    // 7. Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, 17, Eigen::Dynamic, 0, 17, 2* MAX_GNSS_SATELLITES> K_matrix = 
                                state_est_prior.P_matrix * H_matrix.transpose() * 
                                (H_matrix * state_est_prior.P_matrix * H_matrix.transpose() + R_matrix).inverse();

    // 8. Formulate measurement innovations using (14.119)
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 2* MAX_GNSS_SATELLITES, 1> delta_z = Eigen::MatrixXd::Zero(2*gnss_meas.no_meas,1);

    // Range innovations
    delta_z.block(0,0,gnss_meas.no_meas,1) = gnss_meas.meas.block(0,0,gnss_meas.no_meas,1) -
                                                pred_meas.block(0,0,gnss_meas.no_meas,1);
    // Range rates innovations
    delta_z.block(gnss_meas.no_meas,0,gnss_meas.no_meas,1) = gnss_meas.meas.block(0,1,gnss_meas.no_meas,1) -
                                                                pred_meas.block(0,1,gnss_meas.no_meas,1);

    // 9. Update error state estimates using (3.24)
    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,17,1> x_est_new = /*x_est_propagated + */ K_matrix * delta_z;

    // 10. Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,17,17> P_matrix_post = (Eigen::Matrix<double,17,17>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix;

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcefTc state_est_post;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    // Update clock offset + drift estimates
    state_est_post.clock_offset = state_est_prior.clock_offset +  x_est_new(15,0);
    state_est_post.clock_drift = state_est_prior.clock_drift + x_est_new(16,0);

    return state_est_post;
}

StateEstEcefLc lcUpdateKFPosRotEcef (const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcefLc & state_est_prior) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_propagated = Eigen::Matrix<double,15,1>::Zero();
        
    // 5. Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,6,15> H_matrix = Eigen::Matrix<double,6,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position
    H_matrix.block<3,3>(3,0) = - Eigen::Matrix3d::Identity(); // Rotation

    // 6. Set-up measurement noise covariance matrix
    Eigen::Matrix<double, 6,6> R_matrix = pos_rot_meas.cov_mat;

    // 7. Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, 15,6> K_matrix = state_est_prior.P_matrix * H_matrix.transpose() * 
                                (H_matrix * state_est_prior.P_matrix * H_matrix.transpose() + R_matrix).inverse();

    // 8. Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here. See (14.151) for attitude int
    Eigen::Matrix<double,6,1> delta_z;
    delta_z.block<3,1>(0,0) = pos_rot_meas.r_eb_e - state_est_prior.nav_sol.r_eb_e; // pos
    delta_z.block<3,1>(3,0) = deSkew(pos_rot_meas.C_b_e * state_est_prior.nav_sol.C_b_e.transpose() - Eigen::Matrix3d::Identity());// rot

    // 9. Update error state estimates using (3.24)
    Eigen::Matrix<double,15,1> x_est_new = x_est_propagated + K_matrix * delta_z;

    // 10. Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_post = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix;

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcefLc state_est_post;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix = P_matrix_post;
    
    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    return state_est_post;

}

};
