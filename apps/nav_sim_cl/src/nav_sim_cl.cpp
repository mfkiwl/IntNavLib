
#include <iostream>
#include <random>
#include <filesystem>

#include <glog/logging.h>

#include "intnavlib.h"

#include <FGFDMExec.h>

/// @example nav_sim_cl.cpp
/// Integrated navigation demo script - closed loop

using namespace intnavlib;

using Vector3 = Eigen::Matrix<nav_type,3,1>;
using Vector2 = Eigen::Matrix<nav_type,2,1>;
using Matrix3 = Eigen::Matrix<nav_type,3,3>;

static constexpr nav_type kFeetToMeters = 0.3048;

// Helper to get na sol from FDM
NavSolutionNed getNavSolFromFdm(const JSBSim::FGFDMExec & fdm) {

    std::shared_ptr<JSBSim::FGPropagate> prop = fdm.GetPropagate();

    NavSolutionNed true_nav_ned;

    true_nav_ned.time = fdm.GetSimTime();
    true_nav_ned.latitude = prop->GetGeodLatitudeRad();
    true_nav_ned.longitude = prop->GetLongitude();

    // Height: convert to m
    true_nav_ned.height = prop->GetGeodeticAltitude() * kFeetToMeters;

    // Velocity: convert to m/s
    true_nav_ned.v_eb_n = Eigen::Vector3d(prop->GetVel(1), prop->GetVel(2), prop->GetVel(3)) * kFeetToMeters;

    auto C = prop->GetTb2l();
    true_nav_ned.C_b_n << C(1,1), C(1,2), C(1,3),
                            C(2,1), C(2,2), C(2,3),
                            C(3,1), C(3,2), C(3,3);
    
    return true_nav_ned;
}

int main(int argc, char* argv[]) {

    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    if (argc < 2) {
        LOG(ERROR) << "Usage: " << argv[0] << " <path_to_jsbsim_script.xml>\n";
        LOG(ERROR) << "Example: ./build/" << argv[0] << " scripts/c3105.xml\n";
        return 1;
    }

    // Get path from command line
    std::string script_path = argv[1];

    JSBSim::FGFDMExec fdm;

    // Try to load the provided run script
    std::cout << "Loading JSBSim script: " << script_path << std::endl;
    bool ok = fdm.LoadScript(SGPath(script_path));
    if (!ok) {
        LOG(ERROR) << "Failed to load JSBSim script: " << script_path << "\n";
        return 1;
    }

    // Initialize the model according to script
    fdm.RunIC();

    nav_type end_time = 3600; // TODO Get from xml

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

    // Init profile writer
    std::string new_directory = "./results";
    std::string base_filename = std::filesystem::path(script_path).filename().string();
    std::string filename_without_extension = base_filename.substr(0, base_filename.find_last_of('.'));
    std::string eval_filename_out = new_directory + "/" + filename_without_extension + "_eval.csv";
    if (!std::filesystem::exists(new_directory)) {
        std::filesystem::create_directory(new_directory);
    }
    FileWriter eval_data_writer(eval_filename_out);

    // True nav solution
    NavSolutionNed true_nav_ned = getNavSolFromFdm(fdm);
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
    state_est_ecef_init.nav_sol = true_nav_ecef; // Debug: Ideal init
    NavKF nav_filter(state_est_ecef_init, kf_config);

    while (fdm.GetSimTime() < end_time) { // run until script end

        // Run sim step
        fdm.Run();

        // ========= Get ground truth ============

        true_nav_ned = getNavSolFromFdm(fdm);
        true_nav_ecef = nedToEcef(true_nav_ned);

        // LOG(INFO)<< "t=" << true_nav_ned.time
        //           << "  LLA=(" << true_nav_ned.latitude << ", " << true_nav_ned.longitude << ", " << true_nav_ned.height << ")";

        nav_type tor_i = true_nav_ned.time - nav_filter.getTime();

        // ========== IMU Simulation ==========

        // Get true specific force and angular rates
        ImuMeasurements ideal_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // Get imu measurements by applying IMU model
        ImuMeasurements imu_meas = imuModel(ideal_imu_meas, imu_meas_old, imu_errors, tor_i, gen);
        // imu_meas = ideal_imu_meas; // Debug: ideal IMU
        
        // ========== Predict ==========

        nav_filter.lcPredict(imu_meas, tor_i);

        // ========== Update ===========

        nav_type tor_s = true_nav_ned.time - time_last_update;
        if(tor_s >= gnss_config.epoch_interval) {
            PosMeasEcef pos_meas_ecef = genericPosSensModel(true_nav_ecef, 10.0, gen);
            nav_filter.lcUpdatePosEcef(pos_meas_ecef);
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

    LOG(INFO) << "Done!";

    return 0;
}
