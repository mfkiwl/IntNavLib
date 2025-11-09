#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <gtest/gtest.h>
#include <glog/logging.h>

#include "intnavlib.h"

using namespace intnavlib;

using Vector3 = Eigen::Matrix<nav_type,3,1>;
using Vector2 = Eigen::Matrix<nav_type,2,1>;
using Matrix3 = Eigen::Matrix<nav_type,3,3>;

constexpr nav_type max_pos_error = 15.0; // meters

constexpr char test_profile_path[] = "../data/Profile_3.csv";

enum SimType {
    INS_POS,
    INS_POS_ROT,
    INS_GNSS_LC,
    INS_GNSS_TC
};

std::array<SimType, 5> sim_types = {
    SimType::INS_POS,
    SimType::INS_GNSS_LC,
    SimType::INS_GNSS_TC,
    SimType::INS_POS_ROT
};

TEST(navigation_filter, test_integrated)
{ 

    // ============== Init sim ==============

    // Mersenne twister PRNG, initialized with fixed seed for repeatability
    std::mt19937 gen(43);

    // Tactcal grade IMU errors
    ImuErrors imu_errors = tacticalImuErrors();

    // Default GNSS config
    GnssConfig gnss_config = defaultGnssConfig();

    // Tactcal grade IMU - KF config
    KfConfig kf_config = tacticalImuKFConfig();

    // ============ Test each sim type ===========

    for(auto & sim_type : sim_types) {

    // True nav solution
    NavSolutionNed true_nav_ned;
    MotionProfileReader reader(test_profile_path);
    reader.readNextRow(true_nav_ned);
    NavSolutionEcef true_nav_ecef = nedToEcef(true_nav_ned);   
    NavSolutionEcef true_nav_ecef_old = true_nav_ecef;

    // Old IMU measurements
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Vector3::Zero();
    imu_meas_old.quant_residuals_omega = Vector3::Zero();

    // Time of last KF update
    nav_type time_last_update = -1.0;

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned.time,  gnss_config);
    Eigen::Matrix<nav_type, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> gnss_biases = initializeGnssBiases(true_nav_ecef,
                                                                                                        true_nav_ned,
                                                                                                        sat_pos_vel_0,
                                                                                                        gnss_config,
                                                                                                        gen);
    // GNSS measurements at t0
    GnssMeasurements gnss_meas_t0 = generateGnssMeasurements(true_nav_ned.time,
                                                            sat_pos_vel_0,
                                                            true_nav_ned,
                                                            true_nav_ecef,
                                                            gnss_biases, 
                                                            gnss_config,
                                                            gen);

    // Init navigation filter
    StateEstEcef state_est_ecef_init = initStateFromGroundTruth(true_nav_ecef, kf_config, gnss_meas_t0, gen);
    
    NavKF nav_filter(state_est_ecef_init, kf_config);
    
    double pos_error_sum = 0.0;
    long unsigned int count = 0;

    while (reader.readNextRow(true_nav_ned)) {

        // ========= Get ground truth ============

        nav_type tor_i = true_nav_ned.time - nav_filter.getTime();
        true_nav_ecef = nedToEcef(true_nav_ned);

        // ========== IMU Simulation ==========

        // Get true specific force and angular rates
        ImuMeasurements ideal_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);
        // Get imu measurements by applying IMU model
        ImuMeasurements imu_meas = imuModel(ideal_imu_meas, imu_meas_old, imu_errors, tor_i, gen);    

        // ========== Predict ==========

        if (sim_type == SimType::INS_GNSS_TC) {
            nav_filter.tcPredict(imu_meas, tor_i);
        }
        else {
            nav_filter.lcPredict(imu_meas, tor_i);
        }

        // ========== Update =========

        nav_type tor_s = true_nav_ned.time - time_last_update;
        if(tor_s >= gnss_config.epoch_interval) {

            // Simulate GNSS measurements
            SatPosVel sat_pos_vel = satellitePositionsAndVelocities(true_nav_ned.time, gnss_config);
            GnssMeasurements gnss_meas = generateGnssMeasurements(true_nav_ned.time,
                                                                    sat_pos_vel,
                                                                    true_nav_ned,
                                                                    true_nav_ecef,
                                                                    gnss_biases, 
                                                                    gnss_config,
                                                                    gen);
            
            // Loose GNSS Update
            if(sim_type == SimType::INS_GNSS_LC) {
                nav_filter.lcUpdateGnssEcef(gnss_meas, gnss_config);
            }
            // Tight GNSS Update
            else if(sim_type == SimType::INS_GNSS_TC) {
                nav_filter.tcUpdateGnssEcef(gnss_meas, tor_s);
            }
            // Loose position update
            else if(sim_type == SimType::INS_POS) {
                // Simulate position + attitude sensor measurement
                PosMeasEcef pos_meas_ecef = genericPosSensModel(true_nav_ecef, 10.0, gen);
                nav_filter.lcUpdatePosEcef(pos_meas_ecef);
            }
            // Loose position + attitude update
            else if(sim_type == SimType::INS_POS_ROT) {
                // Simulate position + attitude sensor measurement
                PosRotMeasEcef pos_rot_meas_ecef = genericPosRotSensModel(true_nav_ecef, 10.0, 0.01, gen);
                nav_filter.lcUpdatePosRotEcef(pos_rot_meas_ecef);
            }
            time_last_update = true_nav_ned.time;
        }

        // ========== Compute error ==========

        nav_type pos_error = (nav_filter.getStateEst().nav_sol.r_eb_e - true_nav_ecef.r_eb_e).norm();
        pos_error_sum += double(pos_error);
        count++;
        
        // ============= Update simulation state ==============

        true_nav_ecef_old = true_nav_ecef;
        imu_meas_old = imu_meas;
    }

    double avg_position_error = pos_error_sum / double(count);

    EXPECT_LT(avg_position_error, max_pos_error)
        << "Simulation setup n." << sim_type << "failed. \n"
        << "Mean absolute position error should be under "
        << max_pos_error << " m, but was "
        << avg_position_error << " m.";

    // always print summary for visibility
    std::cout << "[INFO] Avg position error (" << sim_type << "): "
            << avg_position_error << " m" << std::endl;

    }

}