#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>

#include "intnavlib.h"

/// @example nav_sim.cpp
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
    std::string errors_sigmas_filename_out = new_directory + "/" + filename_without_extension + "_eval.csv";
    if (!std::filesystem::exists(new_directory)) {
        std::filesystem::create_directory(new_directory);
    }

    // ============== Init sim ==============

    std::random_device rd;
    // Mersenne twister PRNG, initialized with random seed
    std::mt19937 gen(rd());

    // Tactcal grade IMU errors
    ImuErrors imu_errors = tacticalImuErrors();

    // Default GNSS config
    GnssConfig gnss_config = defaultGnssConfig();

    // Tactcal grade IMU - KF config
    KfConfig kf_config = tacticalImuKFConfig();

    // Init profile reader + writers 
    MotionProfileReader reader(motion_profile_filename_in);
    FileWriter eval_data_writer(errors_sigmas_filename_out);

    // True nav solution
    NavSolutionNed true_nav_ned;
    reader.readNextRow(true_nav_ned);
    NavSolutionEcef true_nav_ecef = nedToEcef(true_nav_ned);   
    NavSolutionEcef true_nav_ecef_old = true_nav_ecef;

    // Old IMU measurements
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();

    // Time of last KF update
    double time_last_update = -1.0;

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned.time,  gnss_config);
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> gnss_biases = initializeGnssBiases(true_nav_ecef,
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

    while (reader.readNextRow(true_nav_ned)) {

        // ========= Get ground truth ============

        double tor_i = true_nav_ned.time - nav_filter.getTime();
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

        double tor_s = true_nav_ned.time - time_last_update;
        if(tor_s >= gnss_config.epoch_interval && sim_type != SimType::INS) {

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

        // ========== Write Results ==========

        StateEstEcef state_est_ecef = nav_filter.getStateEst();
        EvalDataEcef eval_data_ecef = getEvalDataEcef(state_est_ecef, true_nav_ecef);
        eval_data_writer.writeEvalDataRow(eval_data_ecef);
        
        // ============= Update simulation state ==============

        true_nav_ecef_old = true_nav_ecef;
        imu_meas_old = imu_meas;
    }

    return 0;
}