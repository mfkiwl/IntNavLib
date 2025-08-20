#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>

#include "intnavlib.h"

/// @example ins_ecef.cpp
/// Integrated navigation demo script

using namespace intnavlib;

int main(int argc, char** argv)
{   
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    if(argc != 3) {
        LOG(ERROR) << "Usage: ./nav_sim <profile_path> <sim_type>";
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

    // Get simulation type
    std::string sim_type(argv[2]);
    if( sim_type != "ins" && 
        sim_type != "ins_pos" && 
        sim_type != "ins_pos_rot" &&
        sim_type != "ins_gnss_lc" && 
        sim_type != "ins_gnss_tc") {
        LOG(ERROR) << "Invalid simulation type: " << sim_type;
        LOG(ERROR) << "Options are: ins, ins_pos_rot, ins_gnss_lc, ins_gnss_tc";
        return 1;
    }

    // ============== Init sim ==============

    std::random_device rd;
    // Mersenne twister PRNG, initialized with random seed
    std::mt19937 gen(rd());

    // Default constructor: tactical grade IMU
    ImuErrors imu_errors;

    GnssConfig gnss_config;

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

    // Current error state estimate
    StateEstEcef state_est_ecef;
    state_est_ecef.valid = true;
    state_est_ecef.nav_sol = true_nav_ecef_old;
    std::normal_distribution att_d{0.0, kf_config.init_att_unc};
    std::normal_distribution vel_d{0.0, kf_config.init_vel_unc};
    std::normal_distribution pos_d{0.0, kf_config.init_pos_unc};
    state_est_ecef.nav_sol.C_b_e = true_nav_ecef_old.C_b_e * rpyToR(Eigen::Vector3d(att_d(gen), att_d(gen), att_d(gen)));
    state_est_ecef.nav_sol.r_eb_e += Eigen::Vector3d(pos_d(gen), pos_d(gen), pos_d(gen));
    state_est_ecef.nav_sol.v_eb_e += Eigen::Vector3d(vel_d(gen), vel_d(gen), vel_d(gen));
    state_est_ecef.acc_bias = Eigen::Vector3d::Zero();
    state_est_ecef.gyro_bias = Eigen::Vector3d::Zero();
    state_est_ecef.clock_offset = 0.0;
    state_est_ecef.clock_drift = 0.0;
    state_est_ecef.P_matrix = initializePMmatrix(kf_config);

    // Time of last KF update
    double time_last_update = 0.0;

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned_t0.time, 
                                                                gnss_config);
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MAX_GNSS_SATELLITES> gnss_biases = initializeGNSSBiases(true_nav_ecef_old,
                                                                                                        true_nav_ned_t0,
                                                                                                        sat_pos_vel_0,
                                                                                                        gnss_config,
                                                                                                        gen);
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

        // ========== Closed loop IMU compensation =========
        
        // correct imu bias using previous state estimation
        imu_meas.f -= state_est_ecef.acc_bias;
        imu_meas.omega -= state_est_ecef.gyro_bias;

        // ========== Predict ==========

        if (sim_type == "ins_gnss_tc") {
            state_est_ecef = tcPredictKF(state_est_ecef, imu_meas, kf_config, tor_i);
        }
        else {
            state_est_ecef = lcPredictKF(state_est_ecef, imu_meas, kf_config, tor_i);
        }

        // ========== Update =========

        double tor_s = true_nav_ned.time - time_last_update;
        if(tor_s >= gnss_config.epoch_interval) {

            time_last_update = true_nav_ned.time;

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
            if(sim_type == "ins_gnss_lc") {
                // Estimate receiver pos + vel with NL LS
                GnssLsPosVelClock pos_vel_clock_gnss_meas_ecef = gnssLsPositionVelocityClock(gnss_meas, state_est_ecef.nav_sol.r_eb_e, state_est_ecef.nav_sol.v_eb_e);
                // Create meas object for loose integration
                GnssPosVelMeasEcef pos_vel_gnss_meas_ecef;
                pos_vel_gnss_meas_ecef.r_ea_e = pos_vel_clock_gnss_meas_ecef.r_ea_e;
                pos_vel_gnss_meas_ecef.v_ea_e = pos_vel_clock_gnss_meas_ecef.v_ea_e;
                pos_vel_gnss_meas_ecef.cov_mat = Eigen::Matrix<double,6,6>::Identity();
                pos_vel_gnss_meas_ecef.cov_mat.block<3,3>(0,0) = pow(gnss_config.lc_pos_sd,2.0) * Eigen::Matrix3d::Identity();
                pos_vel_gnss_meas_ecef.cov_mat.block<3,3>(3,3) = pow(gnss_config.lc_vel_sd,2.0) * Eigen::Matrix3d::Identity();
                // KF update
                state_est_ecef = lcUpdateKFGnssEcef(pos_vel_gnss_meas_ecef, state_est_ecef);
            }
            // Tight GNSS Update
            else if(sim_type == "ins_gnss_tc") {
                state_est_ecef = tcUpdateKFGnssEcef(gnss_meas, state_est_ecef, tor_s);
            }
            // Loose position update
            else if(sim_type == "ins_pos") {
                // Simulate position + attitude sensor measurement
                PosMeasEcef pos_meas_ecef = genericPosSensModel(true_nav_ecef, 10.0, gen);
                state_est_ecef = lcUpdateKFPosEcef(pos_meas_ecef, state_est_ecef);
            }
            // Loose position + attitude update
            else if(sim_type == "ins_pos_rot") {
                // Simulate position + attitude sensor measurement
                PosRotMeasEcef pos_rot_meas_ecef = genericPosRotSensModel(true_nav_ecef, 10.0, 0.01, gen);
                state_est_ecef = lcUpdateKFPosRotEcef(pos_rot_meas_ecef, state_est_ecef);
            }
            // Else no update, pure INS
    
        }

        // ========== Write Results ==========

        ErrorsSigmasEcef errors_sigmas_ecef; 
        errors_sigmas_ecef.time = true_nav_ned.time;
        errors_sigmas_ecef.delta_r_eb_e = state_est_ecef.nav_sol.r_eb_e - true_nav_ecef.r_eb_e;
        errors_sigmas_ecef.delta_v_eb_e = state_est_ecef.nav_sol.v_eb_e - true_nav_ecef.v_eb_e;
        errors_sigmas_ecef.delta_rot_eb_e = deSkew(state_est_ecef.nav_sol.C_b_e * true_nav_ecef.C_b_e.transpose() - Eigen::Matrix3d::Identity());
        
        errors_sigmas_ecef.sigma_delta_r_eb_e << sqrt(state_est_ecef.P_matrix(6,6)), sqrt(state_est_ecef.P_matrix(7,7)), sqrt(state_est_ecef.P_matrix(8,8));
        errors_sigmas_ecef.sigma_delta_v_eb_e << sqrt(state_est_ecef.P_matrix(3,3)), sqrt(state_est_ecef.P_matrix(4,4)), sqrt(state_est_ecef.P_matrix(5,5));
        errors_sigmas_ecef.sigma_delta_rot_eb_e << sqrt(state_est_ecef.P_matrix(0,0)), sqrt(state_est_ecef.P_matrix(1,1)), sqrt(state_est_ecef.P_matrix(2,2));
        errors_sigmas_writer.writeNextRow(errors_sigmas_ecef);
        
        // ============= Update simulation state ==============

        true_nav_ecef_old = true_nav_ecef;
        imu_meas_old = imu_meas;

    }

    return 0;
}