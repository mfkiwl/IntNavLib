#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>

#include "intnavlib.h"

/// @example ins_ecef.cpp
/// Integrated navigation demo script

using namespace intnavlib;

enum SimType {
    INS,
    INS_POS,
    INS_POS_ROT,
    INS_GNSS_LC,
    INS_GNSS_TC,
    UNKNOWN
};

SimType parseSimType(const std::string& sim_type) {
    if (sim_type == "ins") return SimType::INS;
    if (sim_type == "ins_pos")          return SimType::INS_POS;
    if (sim_type == "ins_pos_rot")         return SimType::INS_POS_ROT;
    if (sim_type == "ins_gnss_lc")         return SimType::INS_GNSS_LC;
    if (sim_type == "ins_gnss_tc")         return SimType::INS_GNSS_TC;
    return SimType::UNKNOWN;
}

int main(int argc, char** argv)
{   
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    if(argc != 3) {
        LOG(ERROR) << "Usage: ./nav_sim <profile_path> <sim_type>";
        return 1;
    }

    // Get simulation type
    SimType sim_type = parseSimType(argv[2]);
    if(sim_type == SimType::UNKNOWN) {
        LOG(ERROR) << "Unknown simulation type";
        return 1;
    }

    // Get input profile path and create output dir
    std::string motion_profile_filename_in(argv[1]);
    std::string new_directory = "../results";
    std::string base_filename = std::filesystem::path(motion_profile_filename_in).filename().string();
    std::string filename_without_extension = base_filename.substr(0, base_filename.find_last_of('.'));
    std::string errors_sigmas_filename_out = new_directory + "/" + filename_without_extension + "_errors_sigmas.csv";
    if (!std::filesystem::exists(new_directory)) {
        std::filesystem::create_directory(new_directory);
    }

    // ============== Init sim ==============

    std::random_device rd;
    // Mersenne twister PRNG, initialized with random seed
    std::mt19937 gen(rd());

    // Default constructor: tactical grade IMU
    ImuErrors imu_errors;

    // Default
    GnssConfig gnss_config;

    // Default
    KfConfig kf_config;

    // Init profile reader + writers 
    MotionProfileReader reader(motion_profile_filename_in);
    ErrorsSigmasEcefWriter errors_sigmas_writer(errors_sigmas_filename_out);

    // True old nav solution
    NavSolutionNed true_nav_ned_t0;
    reader.readNextRow(true_nav_ned_t0);
    NavSolutionEcef true_nav_ecef_old = nedToEcef(true_nav_ned_t0);    

    // Old IMU measurements
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();

    // Time of last KF update
    double time_last_update = -1.0;

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned_t0.time,  gnss_config);
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MAX_GNSS_SATELLITES> gnss_biases = initializeGNSSBiases(true_nav_ecef_old,
                                                                                                        true_nav_ned_t0,
                                                                                                        sat_pos_vel_0,
                                                                                                        gnss_config,
                                                                                                        gen);
    // GNSS measurements at t0
    GnssMeasurements gnss_meas_t0 = generateGNSSMeasurements(true_nav_ned_t0.time,
                                                            sat_pos_vel_0,
                                                            true_nav_ned_t0,
                                                            true_nav_ecef_old,
                                                            gnss_biases, 
                                                            gnss_config,
                                                            gen);

    // Init navigation filter state estimate
    StateEstEcef state_est_ecef = initStateFromGroundTruth(true_nav_ecef_old, kf_config, gnss_meas_t0, gen);

    NavSolutionNed true_nav_ned;
    while (reader.readNextRow(true_nav_ned)) {

        // ========= Get ground truth ============

        double tor_i = true_nav_ned.time - state_est_ecef.nav_sol.time;
        NavSolutionEcef true_nav_ecef = nedToEcef(true_nav_ned);

        // ========== IMU Simulation ==========

        // Get true specific force and angular rates
        ImuMeasurements ideal_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);
        // Get imu measurements by applying IMU model
        ImuMeasurements imu_meas = imuModel(ideal_imu_meas, imu_meas_old, imu_errors, tor_i, gen);    

        // ========== Predict ==========

        if (sim_type == SimType::INS_GNSS_TC) {
            state_est_ecef = tcPredictKF(state_est_ecef, imu_meas, kf_config, tor_i);
        }
        else {
            state_est_ecef = lcPredictKF(state_est_ecef, imu_meas, kf_config, tor_i);
        }

        // ========== Update =========

        double tor_s = true_nav_ned.time - time_last_update;
        if(tor_s >= gnss_config.epoch_interval) {

            // Simulate GNSS measurement
            SatPosVel sat_pos_vel = satellitePositionsAndVelocities(true_nav_ned.time, gnss_config);
            GnssMeasurements gnss_meas = generateGNSSMeasurements(true_nav_ned.time,
                                                                    sat_pos_vel,
                                                                    true_nav_ned,
                                                                    true_nav_ecef,
                                                                    gnss_biases, 
                                                                    gnss_config,
                                                                    gen);
            
            // Loose GNSS Update
            if(sim_type == SimType::INS_GNSS_LC) {
                // Estimate receiver pos + vel with NLLS
                GnssPosVelMeasEcef pos_vel_gnss_meas_ecef = gnssLsPositionVelocity(gnss_meas, 
                                                                                    state_est_ecef.nav_sol.r_eb_e, 
                                                                                    state_est_ecef.nav_sol.v_eb_e,
                                                                                    gnss_config);
                state_est_ecef = lcUpdateKFGnssEcef(pos_vel_gnss_meas_ecef, state_est_ecef);
            }
            // Tight GNSS Update
            else if(sim_type == SimType::INS_GNSS_TC) {
                state_est_ecef = tcUpdateKFGnssEcef(gnss_meas, state_est_ecef, tor_s);
            }
            // Loose position update
            else if(sim_type == SimType::INS_POS) {
                // Simulate position + attitude sensor measurement
                PosMeasEcef pos_meas_ecef = genericPosSensModel(true_nav_ecef, 10.0, gen);
                state_est_ecef = lcUpdateKFPosEcef(pos_meas_ecef, state_est_ecef);
            }
            // Loose position + attitude update
            else if(sim_type == SimType::INS_POS_ROT) {
                // Simulate position + attitude sensor measurement
                PosRotMeasEcef pos_rot_meas_ecef = genericPosRotSensModel(true_nav_ecef, 10.0, 0.01, gen);
                state_est_ecef = lcUpdateKFPosRotEcef(pos_rot_meas_ecef, state_est_ecef);
            }
            // Else no update, pure INS
            else continue;

            time_last_update = true_nav_ned.time;
        }

        // ========== Write Results ==========

        ErrorsSigmasEcef errors_sigmas_ecef = getErrorsSigmasEcef(state_est_ecef, true_nav_ecef);
        errors_sigmas_writer.writeNextRow(errors_sigmas_ecef);
        
        // ============= Update simulation state ==============

        true_nav_ecef_old = true_nav_ecef;
        imu_meas_old = imu_meas;

    }

    return 0;
}