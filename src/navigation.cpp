#include "navigation.h"

namespace intnavlib {

NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i) {
    NavSolutionEcef new_nav;
    new_nav.time = imu_meas.time;
    // ATTITUDE UPDATE
    // From (2.145) determine the Earth rotation over the update interval
    // C_Earth = C_e_i' * old_C_e_i
    double alpha_ie = kOmega_ie * tor_i;
    double cos_alpha_ie = cos(alpha_ie);
    double sin_alpha_ie = sin(alpha_ie);
    Eigen::Matrix3d C_Earth;
    C_Earth << cos_alpha_ie, sin_alpha_ie, 0.0,
                -sin_alpha_ie, cos_alpha_ie, 0.0,
                            0.0,             0.0,  1.0;
    // Calculate attitude increment, magnitude, and skew-symmetric matrix
    Eigen::Vector3d alpha_ib_b = imu_meas.omega * tor_i;
    double mag_alpha = alpha_ib_b.norm();
    Eigen::Matrix3d Alpha_ib_b = skewSymmetric(alpha_ib_b);  
    // Obtain coordinate transformation matrix from the new attitude w.r.t. an
    // inertial frame to the old using Rodrigues' formula, (5.73)
    Eigen::Matrix3d C_new_old;
    if (mag_alpha>1.0e-8) {
        C_new_old = Eigen::Matrix3d::Identity() + 
                    ((sin(mag_alpha) / mag_alpha) * Alpha_ib_b) +
                    ((1.0 - cos(mag_alpha)) / pow(mag_alpha,2.0)) * Alpha_ib_b * Alpha_ib_b;
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
    // Transform specific force to ECEF-frame resolving axes using (5.85)
    Eigen::Vector3d f_ib_e = ave_C_b_e * imu_meas.f;
    // UPDATE VELOCITY
    // From (5.36)
    Eigen::Vector3d kOmega_ie_vec;
    kOmega_ie_vec << 0.0 , 0.0 , kOmega_ie;   
    new_nav.v_eb_e = old_nav.v_eb_e + tor_i * (f_ib_e + gravityEcef(old_nav.r_eb_e) -
                            2.0 * skewSymmetric(kOmega_ie_vec) * old_nav.v_eb_e);
    // UPDATE CARTESIAN POSITION
    // From (5.38),
    new_nav.r_eb_e = old_nav.r_eb_e + (new_nav.v_eb_e + old_nav.v_eb_e) * 0.5 * tor_i; 
    return new_nav;
}

StateEstEcef lcPredictKF(const StateEstEcef & state_est_old, 
                        const ImuMeasurements & imu_meas,
                        const KfConfig & lc_kf_config,
                        const double & tor_i) {
    // Compensate IMU measurements
    ImuMeasurements imu_meas_comp = imu_meas;
    imu_meas_comp.f -= state_est_old.acc_bias;
    imu_meas_comp.omega -= state_est_old.gyro_bias;
    // Prop uncertainty
    StateEstEcef state_est_ecef = state_est_old;
    state_est_ecef.P_matrix.block<15,15>(0,0) = lcPropUnc(state_est_old.P_matrix.block<15,15>(0,0), 
                                                            state_est_old.nav_sol,
                                                            ecefToNed(state_est_old.nav_sol),
                                                            imu_meas_comp,
                                                            lc_kf_config,
                                                            tor_i);
    // Predict state
    state_est_ecef.nav_sol = navEquationsEcef(state_est_ecef.nav_sol, imu_meas_comp, tor_i);
    return state_est_ecef;
}

StateEstEcef tcPredictKF(const StateEstEcef & state_est_old, 
                        const ImuMeasurements & imu_meas,
                        const KfConfig & kf_config,
                        const double & tor_i) {
    // Compensate IMU measurements
    ImuMeasurements imu_meas_comp = imu_meas;
    imu_meas_comp.f -= state_est_old.acc_bias;
    imu_meas_comp.omega -= state_est_old.gyro_bias;
    // Prop uncertainty
    StateEstEcef state_est_ecef = state_est_old;
    state_est_ecef.P_matrix = tcPropUnc(state_est_old.P_matrix, 
                                        state_est_old.nav_sol,
                                        ecefToNed(state_est_old.nav_sol),
                                        imu_meas_comp,
                                        kf_config,
                                        tor_i);
    // Predict state
    state_est_ecef.nav_sol = navEquationsEcef(state_est_ecef.nav_sol, imu_meas_comp, tor_i);
    return state_est_ecef;
}

Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & lc_kf_config,
                                        const double & tor_i) {

    double geocentric_radius = kR0 / sqrt(1.0 - pow(kEccentricity * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - kEccentricity*kEccentricity, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Eigen::Vector3d kOmega_ie_vec;
    kOmega_ie_vec << 0.0,0.0,kOmega_ie;
    Eigen::Matrix3d Omega_ie = skewSymmetric(kOmega_ie_vec);

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
    Q_prime_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * lc_kf_config.gyro_noise_psd * tor_i;
    Q_prime_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * lc_kf_config.accel_noise_psd * tor_i;
    Q_prime_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * lc_kf_config.accel_bias_psd * tor_i;
    Q_prime_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * lc_kf_config.gyro_bias_psd * tor_i;

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
    double geocentric_radius = kR0 / sqrt(1.0 - pow(kEccentricity * sin(old_nav_est_ned.latitude),2.0)) *
        sqrt(pow(cos(old_nav_est_ned.latitude), 2.0) + pow(1.0 - kEccentricity*kEccentricity, 2.0) * pow(sin(old_nav_est_ned.latitude), 2.0)); // from (2.137)

    // Skew symmetric matrix of Earth rate
    Eigen::Vector3d kOmega_ie_vec;
    kOmega_ie_vec << 0.0,0.0,kOmega_ie;
    Eigen::Matrix3d Omega_ie = skewSymmetric(kOmega_ie_vec);

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
    Q_prime_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * tc_kf_config.gyro_noise_psd * tor_i;
    Q_prime_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * tc_kf_config.accel_noise_psd * tor_i;
    Q_prime_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * tc_kf_config.accel_bias_psd * tor_i;
    Q_prime_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * tc_kf_config.gyro_bias_psd * tor_i;

    // Clock offset, drift
    Q_prime_matrix(15,15) = tc_kf_config.clock_phase_psd * tor_i;
    Q_prime_matrix(16,16) = tc_kf_config.clock_freq_psd * tor_i;

    // 4. Propagate state estimation error covariance matrix using (3.46)
    Eigen::Matrix<double,17,17> P_matrix = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
    return P_matrix;
}

StateEstEcef lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                const StateEstEcef & state_est_prior,
                                const double & p_value) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_propagated = Eigen::Matrix<double,15,1>::Zero();
        
    // Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,3,15> H_matrix = Eigen::Matrix<double,3,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position

    // Set-up measurement noise covariance matrix
    Eigen::Matrix3d R_matrix = pos_meas.cov_mat;

    // Calculate Kalman gain using (3.21)
    Eigen::Matrix3d S_matrix_inv = (H_matrix * state_est_prior.P_matrix.block<15,15>(0,0) * H_matrix.transpose() + R_matrix).inverse();
    Eigen::Matrix<double, 15,3> K_matrix = state_est_prior.P_matrix.block<15,15>(0,0) * H_matrix.transpose() * S_matrix_inv;

    // Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Eigen::Vector3d delta_z = pos_meas.r_eb_e - state_est_prior.nav_sol.r_eb_e;

    // Update error state estimates using (3.24)
    Eigen::Matrix<double,15,1> x_est_new = x_est_propagated + K_matrix * delta_z;

    // Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_post = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix.block<15,15>(0,0);

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix.block<15,15>(0,0) = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    // Real-time consistency check
    // See: Estimation with Applications to Tracking and Navigation -- Yaakov Bar-Shalom et al. p.237
    // The quadratic form res' * inv(S) * res is a determination of the standard Chi squared distribution with one dof.
    // We check it's inside the desired p_value
    // Simplified version: no history, just current measurement, 1 * n_res dof's.
    double chi2 = delta_z.dot(S_matrix_inv * delta_z);
    unsigned int n_dofs = delta_z.rows();
    boost::math::chi_squared chi2_dist(n_dofs);
    double chi2_thresh = boost::math::quantile(chi2_dist, p_value);
    state_est_post.valid = true;
    if (chi2 > chi2_thresh)
        // Set invalid update flag
        state_est_post.valid = false;
    
    // Return 
    return state_est_post;
}

StateEstEcef lcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const GnssConfig & gnss_config,
                                    const double & p_value) {

    // Get position + velocity estimate with nlls
    GnssPosVelMeasEcef pos_vel_gnss_meas = gnssLsPositionVelocity(gnss_meas, 
                                                                state_est_prior.nav_sol.r_eb_e, 
                                                                state_est_prior.nav_sol.v_eb_e,
                                                                gnss_config);
        
    // Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,6,15> H_matrix = Eigen::Matrix<double,6,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position
    H_matrix.block<3,3>(3,3) = - Eigen::Matrix3d::Identity(); // Velocity

    // Set-up measurement noise covariance matrix
    Eigen::Matrix<double,6,6> R_matrix = pos_vel_gnss_meas.cov_mat;

    // Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, 6,6> S_matrix_inv = (H_matrix * state_est_prior.P_matrix.block<15,15>(0,0) * H_matrix.transpose() + R_matrix).inverse();
    Eigen::Matrix<double, 15,6> K_matrix = state_est_prior.P_matrix.block<15,15>(0,0) * H_matrix.transpose() * S_matrix_inv;

    // Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here
    Eigen::Matrix<double,6,1> delta_z;
    delta_z.block<3,1>(0,0) = pos_vel_gnss_meas.r_ea_e - state_est_prior.nav_sol.r_eb_e;
    delta_z.block<3,1>(3,0) = pos_vel_gnss_meas.v_ea_e - state_est_prior.nav_sol.v_eb_e;

    // Update error state estimates using (3.24)
    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_new = /*x_est_propagated + */ K_matrix * delta_z;

    // Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_post = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix.block<15,15>(0,0);

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix.block<15,15>(0,0) = P_matrix_post;

    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    // Real-time consistency check
    // See: Estimation with Applications to Tracking and Navigation -- Yaakov Bar-Shalom et al. p.237
    // The quadratic form res' * inv(S) * res is a determination of the standard Chi squared distribution with one dof.
    // We check it's inside the desired p_value
    // Simplified version: no history, just current measurement, 1 * n_res dof's.
    double chi2 = delta_z.dot(S_matrix_inv * delta_z);
    unsigned int n_dofs = delta_z.rows();
    boost::math::chi_squared chi2_dist(n_dofs);
    double chi2_thresh = boost::math::quantile(chi2_dist, p_value);
    state_est_post.valid = true;
    if (chi2 > chi2_thresh)
        // Set invalid update flag
        state_est_post.valid = false;

    return state_est_post;
}


StateEstEcef tcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const double & tor_s,
                                    const double & p_value) {

    // Compute predicted range and range rate measurements from prior state est
    Eigen::Matrix<double, Eigen::Dynamic, 2, 0, kMaxGnssSatellites, 2> pred_meas(gnss_meas.no_meas,2);
    // Line of sight unit vectors in ecef
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, kMaxGnssSatellites, 3> u_as_e_T(gnss_meas.no_meas,3);

    // For each satellite
    for(int j=0; j<gnss_meas.no_meas; j++) {

        // Predicted range, approximated
        Eigen::Vector3d delta_r = gnss_meas.meas.block<1,3>(j,2).transpose() - state_est_prior.nav_sol.r_eb_e;
        double range = delta_r.norm();

        // Calculate ecef rotation during signal transit time using (8.36)
        Eigen::Matrix3d C_e_I;
        C_e_I << 1, kOmega_ie * range / c, 0,
                -kOmega_ie * range / c, 1, 0,
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
        Eigen::Vector3d kOmega_ie_vec;
        kOmega_ie_vec << 0.0,0.0,kOmega_ie;
        Eigen::Matrix3d Omega_ie = skewSymmetric(kOmega_ie_vec);
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
    Eigen::Matrix<double, Eigen::Dynamic, 17, 0, 2 * kMaxGnssSatellites, 17> H_matrix = 
        Eigen::Matrix<double, Eigen::Dynamic, 17, 0, 2 * kMaxGnssSatellites, 17>::Zero(2*gnss_meas.no_meas, 17);

    // Ranges
    H_matrix.block(0,6,gnss_meas.no_meas,3) = u_as_e_T.block(0,0,gnss_meas.no_meas,3);
    H_matrix.block(0,15,gnss_meas.no_meas,1) = Eigen::MatrixXd::Ones(gnss_meas.no_meas,1);
    // Range rates
    H_matrix.block(gnss_meas.no_meas,3, gnss_meas.no_meas,3) = u_as_e_T.block(0,0,gnss_meas.no_meas, 3);
    H_matrix.block(gnss_meas.no_meas,16, gnss_meas.no_meas,1) = Eigen::MatrixXd::Ones(gnss_meas.no_meas,1);

    
    // 6. Set-up measurement noise covariance matrix
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites> R_matrix = gnss_meas.cov_mat;

    // 7. Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites> S_matrix_inv = (H_matrix * state_est_prior.P_matrix * H_matrix.transpose() + R_matrix).inverse();
    Eigen::Matrix<double, 17, Eigen::Dynamic, 0, 17, 2* kMaxGnssSatellites> K_matrix = state_est_prior.P_matrix * H_matrix.transpose() * S_matrix_inv;

    // 8. Formulate measurement innovations using (14.119)
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, 2* kMaxGnssSatellites, 1> delta_z = Eigen::MatrixXd::Zero(2*gnss_meas.no_meas,1);

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
    StateEstEcef state_est_post;
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

    // Real-time consistency check
    // See: Estimation with Applications to Tracking and Navigation -- Yaakov Bar-Shalom et al. p.237
    // The quadratic form res' * inv(S) * res is a determination of the standard Chi squared distribution with one dof.
    // We check it's inside the desired p_value
    // Simplified version: no history, just current measurement, 1 * n_res dof's.
    double chi2 = delta_z.dot(S_matrix_inv * delta_z);
    unsigned int n_dofs = delta_z.rows();
    boost::math::chi_squared chi2_dist(n_dofs);
    double chi2_thresh = boost::math::quantile(chi2_dist, p_value);
    state_est_post.valid = true;
    if (chi2 > chi2_thresh)
        // Set invalid update flag
        state_est_post.valid = false;

    return state_est_post;
}

StateEstEcef lcUpdateKFPosRotEcef (const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const double & p_value) {

    // A priori error state is always zero in closed loop filter
    Eigen::Matrix<double,15,1> x_est_propagated = Eigen::Matrix<double,15,1>::Zero();
        
    // 5. Set-up measurement matrix using (14.115)
    Eigen::Matrix<double,6,15> H_matrix = Eigen::Matrix<double,6,15>::Zero();
    H_matrix.block<3,3>(0,6) = - Eigen::Matrix3d::Identity(); // Position
    H_matrix.block<3,3>(3,0) = - Eigen::Matrix3d::Identity(); // Rotation

    // 6. Set-up measurement noise covariance matrix
    Eigen::Matrix<double, 6,6> R_matrix = pos_rot_meas.cov_mat;

    // 7. Calculate Kalman gain using (3.21)
    Eigen::Matrix<double, 6,6> S_matrix_inv = (H_matrix * state_est_prior.P_matrix.block<15,15>(0,0) * H_matrix.transpose() + R_matrix).inverse();
    Eigen::Matrix<double, 15,6> K_matrix = state_est_prior.P_matrix.block<15,15>(0,0) * H_matrix.transpose() * S_matrix_inv;

    // 8. Formulate measurement innovations using (14.102), noting that zero
    // lever arm is assumed here. See (14.151) for attitude int
    Eigen::Matrix<double,6,1> delta_z;
    delta_z.block<3,1>(0,0) = pos_rot_meas.r_eb_e - state_est_prior.nav_sol.r_eb_e; // pos
    delta_z.block<3,1>(3,0) = deSkew(pos_rot_meas.C_b_e * state_est_prior.nav_sol.C_b_e.transpose() - Eigen::Matrix3d::Identity());// rot

    // 9. Update error state estimates using (3.24)
    Eigen::Matrix<double,15,1> x_est_new = x_est_propagated + K_matrix * delta_z;

    // 10. Update state estimation error covariance matrix using (3.25)
    Eigen::Matrix<double,15,15> P_matrix_post = (Eigen::Matrix<double,15,15>::Identity() 
                                            - K_matrix * H_matrix) * state_est_prior.P_matrix.block<15,15>(0,0);

    // CLOSED-LOOP CORRECTION

    // Correct attitude, velocity, and position using (14.7-9)
    StateEstEcef state_est_post;
    state_est_post.nav_sol.time = state_est_prior.nav_sol.time;
    state_est_post.nav_sol.C_b_e = (Eigen::Matrix3d::Identity() - 
                                skewSymmetric(x_est_new.block<3,1>(0,0))) * 
                                state_est_prior.nav_sol.C_b_e;
    state_est_post.nav_sol.v_eb_e = state_est_prior.nav_sol.v_eb_e - x_est_new.block<3,1>(3,0);
    state_est_post.nav_sol.r_eb_e = state_est_prior.nav_sol.r_eb_e - x_est_new.block<3,1>(6,0);
    state_est_post.P_matrix.block<15,15>(0,0) = P_matrix_post;
    
    // Update IMU bias estimates
    state_est_post.acc_bias = state_est_prior.acc_bias + x_est_new.block<3,1>(9,0);
    state_est_post.gyro_bias = state_est_prior.gyro_bias + x_est_new.block<3,1>(12,0);

    // Real-time consistency check
    // See: Estimation with Applications to Tracking and Navigation -- Yaakov Bar-Shalom et al. p.237
    // The quadratic form res' * inv(S) * res is a determination of the standard Chi squared distribution with one dof.
    // We check it's inside the desired p_value
    // Simplified version: no history, just current measurement, 1 * n_res dof's.
    double chi2 = delta_z.dot(S_matrix_inv * delta_z);
    unsigned int n_dofs = delta_z.rows();
    boost::math::chi_squared chi2_dist(n_dofs);
    double chi2_thresh = boost::math::quantile(chi2_dist, p_value);
    state_est_post.valid = true;
    if (chi2 > chi2_thresh)
        // Set invalid update flag
        state_est_post.valid = false;

    return state_est_post;

}

GnssLsPosVelClock gnssLsPositionVelocityClock(const GnssMeasurements & gnss_measurements,
                                            const Eigen::Vector3d & prior_r_ea_e,
                                            const Eigen::Vector3d & prior_v_ea_e) {
    GnssLsPosVelClock est_pos_vel;
    // POSITION AND CLOCK OFFSET
    Eigen::Vector4d x_pred, x_est;
    x_pred.segment<3>(0) = prior_r_ea_e;
    x_pred(3) = 0;
    double test_convergence = 1.0;
    while (test_convergence > 0.0001) {
        Eigen::Matrix<double, Eigen::Dynamic, 4, 0, kMaxGnssSatellites> H_matrix(gnss_measurements.no_meas, 4);
        Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> pred_meas(gnss_measurements.no_meas);
        for (int j = 0; j < gnss_measurements.no_meas; ++j) {
            Eigen::Vector3d delta_r = gnss_measurements.meas.block<1,3>(j,2).transpose() - x_pred.segment<3>(0);
            double approx_range = delta_r.norm();
            Eigen::Matrix3d C_e_I;
            C_e_I << 1, kOmega_ie * approx_range / c, 0,
                     -kOmega_ie * approx_range / c, 1, 0,
                     0, 0, 1;
            delta_r = C_e_I * gnss_measurements.meas.block<1,3>(j,2).transpose() - x_pred.segment<3>(0);
            double range = delta_r.norm();
            pred_meas(j) = range + x_pred(3);
            H_matrix.block<1,3>(j,0) = -delta_r.transpose() / range;
            H_matrix(j,3) = 1;
        }
        x_est = x_pred + (H_matrix.transpose() * H_matrix).inverse() * H_matrix.transpose() * (gnss_measurements.meas.col(0).head(gnss_measurements.no_meas) - pred_meas.head(gnss_measurements.no_meas));
        test_convergence = (x_est - x_pred).norm();
        x_pred = x_est;
    }

    est_pos_vel.r_ea_e = x_est.segment<3>(0);
    est_pos_vel.clock(0) = x_est(3);
    // VELOCITY AND CLOCK DRIFT
    Eigen::Matrix3d Omega_ie = skewSymmetric(Eigen::Vector3d(0, 0, kOmega_ie));
    x_pred.segment<3>(0) = prior_v_ea_e;
    x_pred(3) = 0;
    test_convergence = 1.0;
    while (test_convergence > 0.0001) {
        Eigen::Matrix<double, Eigen::Dynamic, 4, 0, kMaxGnssSatellites> H_matrix(gnss_measurements.no_meas, 4);
        Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> pred_meas(gnss_measurements.no_meas);
        for (int j = 0; j < gnss_measurements.no_meas; ++j) {
            Eigen::Vector3d delta_r = gnss_measurements.meas.block<1,3>(j,2).transpose() - est_pos_vel.r_ea_e;
            double approx_range = delta_r.norm();
            Eigen::Matrix3d C_e_I;
            C_e_I << 1, kOmega_ie * approx_range / c, 0,
                     -kOmega_ie * approx_range / c, 1, 0,
                     0, 0, 1;
            delta_r = C_e_I * gnss_measurements.meas.block<1,3>(j,2).transpose() - est_pos_vel.r_ea_e;
            double range = delta_r.norm();
            Eigen::Vector3d u_as_e = delta_r / range;
            Eigen::Vector3d sat_velocity = gnss_measurements.meas.block<1,3>(j,5).transpose();
            Eigen::Vector3d sat_position = gnss_measurements.meas.block<1,3>(j,2).transpose();
            double range_rate = u_as_e.transpose() * (C_e_I * (sat_velocity + Omega_ie * sat_position) - (x_pred.segment<3>(0) + Omega_ie * est_pos_vel.r_ea_e));
            pred_meas(j) = range_rate + x_pred(3);
            H_matrix.block<1,3>(j,0) = -u_as_e.transpose();
            H_matrix(j,3) = 1;
        }
        x_est = x_pred + (H_matrix.transpose() * H_matrix).inverse() * H_matrix.transpose() * (gnss_measurements.meas.col(1).head(gnss_measurements.no_meas) - pred_meas.head(gnss_measurements.no_meas));
        test_convergence = (x_est - x_pred).norm();
        x_pred = x_est;
    }
    est_pos_vel.v_ea_e = x_est.segment<3>(0);
    est_pos_vel.clock(1) = x_est(3);
    return est_pos_vel;
}

GnssPosVelMeasEcef gnssLsPositionVelocity(const GnssMeasurements & gnss_measurements,
                                            const Eigen::Vector3d & prior_r_ea_e,
                                            const Eigen::Vector3d & prior_v_ea_e,
                                            const GnssConfig & gnss_config) {
    
    GnssLsPosVelClock est_pos_vel = gnssLsPositionVelocityClock(gnss_measurements, prior_r_ea_e, prior_v_ea_e);
    GnssPosVelMeasEcef pos_vel_gnss_meas_ecef;
    pos_vel_gnss_meas_ecef.r_ea_e = est_pos_vel.r_ea_e;
    pos_vel_gnss_meas_ecef.v_ea_e = est_pos_vel.v_ea_e;
    pos_vel_gnss_meas_ecef.cov_mat = Eigen::Matrix<double,6,6>::Identity();
    pos_vel_gnss_meas_ecef.cov_mat.block<3,3>(0,0) = pow(gnss_config.lc_pos_sd,2.0) * Eigen::Matrix3d::Identity();
    pos_vel_gnss_meas_ecef.cov_mat.block<3,3>(3,3) = pow(gnss_config.lc_vel_sd,2.0) * Eigen::Matrix3d::Identity();
    return pos_vel_gnss_meas_ecef;
}

Eigen::Matrix<double,17,17> initializePMmatrix(const KfConfig & kf_config) {
    Eigen::Matrix<double,17,17> P_matrix;
    // Initialize error covariance matrix
    P_matrix = Eigen::Matrix<double,17,17>::Zero();
    P_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * pow(kf_config.init_att_unc,2); // attitude error
    P_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * pow(kf_config.init_vel_unc,2); // vel error
    P_matrix.block<3,3>(6,6) = Eigen::Matrix3d::Identity() * pow(kf_config.init_pos_unc,2); // pos error
    P_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * pow(kf_config.init_b_a_unc,2); // acc bias error
    P_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * pow(kf_config.init_b_g_unc,2); // gyro bias error
    P_matrix(15,15) = pow(kf_config.init_clock_offset_unc,2); // clock offset error
    P_matrix(16,16) = pow(kf_config.init_clock_drift_unc,2); // clock drift error
    return P_matrix;
}

StateEstEcef initStateFromGroundTruth(const NavSolutionEcef & true_nav_ecef, const KfConfig & kf_config, const GnssMeasurements & gnss_meas, std::mt19937 & gen) {
    StateEstEcef state_est_ecef;
    state_est_ecef.valid = true;
    state_est_ecef.nav_sol = true_nav_ecef;
    std::normal_distribution att_d{0.0, kf_config.init_att_unc};
    std::normal_distribution vel_d{0.0, kf_config.init_vel_unc};
    std::normal_distribution pos_d{0.0, kf_config.init_pos_unc};
    state_est_ecef.nav_sol.C_b_e = true_nav_ecef.C_b_e * rpyToR(Eigen::Vector3d(att_d(gen), att_d(gen), att_d(gen)));
    state_est_ecef.nav_sol.r_eb_e += Eigen::Vector3d(pos_d(gen), pos_d(gen), pos_d(gen));
    state_est_ecef.nav_sol.v_eb_e += Eigen::Vector3d(vel_d(gen), vel_d(gen), vel_d(gen));
    state_est_ecef.acc_bias = Eigen::Vector3d::Zero();
    state_est_ecef.gyro_bias = Eigen::Vector3d::Zero();
    state_est_ecef.P_matrix = initializePMmatrix(kf_config);
    // Get clock states using NLLS solver
    GnssLsPosVelClock gnss_pos_vel_clock_est = gnssLsPositionVelocityClock(gnss_meas, state_est_ecef.nav_sol.r_eb_e, state_est_ecef.nav_sol.v_eb_e);
    state_est_ecef.clock_offset = gnss_pos_vel_clock_est.clock(0);
    state_est_ecef.clock_drift = gnss_pos_vel_clock_est.clock(1);
    return state_est_ecef;
}

};
