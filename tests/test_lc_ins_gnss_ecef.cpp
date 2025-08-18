#include <gtest/gtest.h>
#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>

#include "intnavlib.h"

using namespace intnavlib;

constexpr double max_pos_error = 15.0; // meters
constexpr char test_profile_path[] = "../data/Profile_3.csv";  // Provide your test profile path

TEST(ins_gnss, test_tc_ins_gnss_ecef)
{   
    
    // Init motion profile reader & writer
    MotionProfileReader reader(test_profile_path);

    // ============== Random gen ==============

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(42);

    // ============== IMU parameters ==============

    ImuErrors imu_errors;

    // ============== GNSS config ==============

    GnssConfig gnss_config;

    // ============== KF config ==============

    KfConfig lc_kf_config;

    // ============ Declare persistent variables ============

    // Navigation solutions
    // Ground truth nav solution in ned
    NavSolutionNed true_nav_ned;
    NavSolutionNed true_nav_ned_old;
    // Ground truth nav solution in ecef
    NavSolutionEcef true_nav_ecef;
    NavSolutionEcef true_nav_ecef_old;
    // Estimated nav solution in ecef
    NavSolutionEcef est_nav_ecef;
    // Estimated nav solution in ned
    NavSolutionNed est_nav_ned;
    // Estimated biases
    Eigen::Vector3d est_acc_bias = Eigen::Vector3d::Zero();
    Eigen::Vector3d est_gyro_bias = Eigen::Vector3d::Zero();

    // IMU measurements
    // Ground truth imu measurements from kinematics
    ImuMeasurements true_imu_meas;
    // Simulated imu measurements
    ImuMeasurements imu_meas;
    // Old simulated imu measurements
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();

    // Pos sensor measurements
    PosMeasEcef pos_meas_ecef;

    // Error state uncertainty
    Eigen::Matrix<double,15,15> P_matrix = initializeLcPMmatrix(lc_kf_config);

    // Init nav solution
    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);

    est_nav_ecef = true_nav_ecef_old;
    est_nav_ecef.r_eb_e += 10 * Eigen::Vector3d::Ones();
    est_nav_ecef.v_eb_e += 0.1 * Eigen::Vector3d::Ones();

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned_old.time, 
                                                            gnss_config);
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MAX_GNSS_SATELLITES> gnss_biases = 
                                    initializeGNSSBiases(true_nav_ecef_old,
                                                        true_nav_ned_old,
                                                        sat_pos_vel_0,
                                                        gnss_config,
                                                        gen);
    
    // Times
    // Current time - last time
    double tor_i;
    double time_last_gnss = 0.0;
    
    double pos_error_sum = 0.0;
    long unsigned int count = 0;
    while (reader.readNextRow(true_nav_ned)) {

        tor_i = true_nav_ned.time - true_nav_ned_old.time;
        
        // Get true ecef nav sol from motion profile in ned
        true_nav_ecef = nedToEcef(true_nav_ned);

        // ========== IMU SIMULATION ==========

        // Get true specific force and angular rates
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);
        // Get imu measurements by applying IMU model
        imu_meas = imuModel(true_imu_meas, imu_meas_old, imu_errors, tor_i, gen);
        // correct imu bias using previous state estimation
        imu_meas.f -= est_acc_bias;
        imu_meas.omega -= est_gyro_bias;

        // ========== NAV EQUATIONS ==========

        // Predict ecef nav solution (INS)
        est_nav_ecef = navEquationsEcef(est_nav_ecef, imu_meas, tor_i);
        est_nav_ned = ecefToNed(est_nav_ecef);

        // ========== PROP UNCERTAINTIES ==========
        
        // Groves actually puts this in the update part, with a large dt. 
        // But the the underlying approximation hold only
        // if dt is low enough. Therefore slower but better to put it in the prop stage.
        P_matrix  = lcPropUnc(P_matrix, 
                            est_nav_ecef,
                            est_nav_ned,
                            imu_meas,
                            lc_kf_config,
                            tor_i);

        // ========== INTEGRATE  MEASUREMENTS ==========

        
        double tor_s = true_nav_ned.time - time_last_gnss;
        if( tor_s >= gnss_config.epoch_interval) {

            time_last_gnss = true_nav_ned.time;

            // Simulate GNSS measurement
            SatPosVel sat_pos_vel = satellitePositionsAndVelocities(true_nav_ned.time, gnss_config);
            GnssMeasurements gnss_meas = generateGNSSMeasurements(true_nav_ned.time,
                                                                sat_pos_vel,
                                                                true_nav_ned,
                                                                true_nav_ecef,
                                                                gnss_biases, 
                                                                gnss_config,
                                                                gen);

            // Estimate receiver pos + vel with NL LS
            GnssLsPosVelClock pos_vel_clock_gnss_meas_ecef = gnssLsPositionVelocityClock(gnss_meas,
                                                                        est_nav_ecef.r_eb_e,
                                                                        est_nav_ecef.v_eb_e);

            // Create meas object for integration
            GnssPosVelMeasEcef pos_vel_gnss_meas_ecef;
            pos_vel_gnss_meas_ecef.r_ea_e = pos_vel_clock_gnss_meas_ecef.r_ea_e;
            pos_vel_gnss_meas_ecef.v_ea_e = pos_vel_clock_gnss_meas_ecef.v_ea_e;
            pos_vel_gnss_meas_ecef.cov_mat = Eigen::Matrix<double,6,6>::Identity();
            pos_vel_gnss_meas_ecef.cov_mat.block<3,3>(0,0) = pow(gnss_config.lc_pos_sd,2.0) * Eigen::Matrix3d::Identity();
            pos_vel_gnss_meas_ecef.cov_mat.block<3,3>(3,3) = pow(gnss_config.lc_vel_sd,2.0) * Eigen::Matrix3d::Identity();


            // KF update -> update posterior
            // if no update, best est is prior
            StateEstEcefLc est_state_ecef_prior;

            est_state_ecef_prior.P_matrix = P_matrix;
            est_state_ecef_prior.nav_sol = est_nav_ecef;
            est_state_ecef_prior.acc_bias = est_acc_bias;
            est_state_ecef_prior.gyro_bias = est_gyro_bias;

            StateEstEcefLc est_state_ecef_post = lcUpdateKFGnssEcef(pos_vel_gnss_meas_ecef, 
                                                                    est_state_ecef_prior);

            P_matrix = est_state_ecef_post.P_matrix;
            est_nav_ecef = est_state_ecef_post.nav_sol;
            est_acc_bias = est_state_ecef_post.acc_bias;
            est_gyro_bias = est_state_ecef_post.gyro_bias;

        }
        
        // ========== COMPUTE ERRORS ==========

        double pos_error = (est_nav_ecef.r_eb_e - true_nav_ecef.r_eb_e).norm();
        pos_error_sum += pos_error;
        count++;

        true_nav_ecef_old = true_nav_ecef;
        true_nav_ned_old = true_nav_ned;
        imu_meas_old = imu_meas;
    }

    double avg_position_error = pos_error_sum / (double) count;

    EXPECT_LT(avg_position_error, max_pos_error)
        << "Mean absolute position error should be under "
        << max_pos_error << " m, but was "
        << avg_position_error << " m.";

}