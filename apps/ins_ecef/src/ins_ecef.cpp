#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>


#include "intnavlib.h"

/// @example ins_ecef.cpp
/// Inertial navigation example in ECEF frame

using namespace intnavlib;

int main(int argc, char** argv)
{   
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    if(argc != 2) {
        LOG(ERROR) << "Pass the motion profile path";
        return 1;
    }

    // Input profile filename
    std::string motion_profile_filename_in(argv[1]);
    
    // Output profile filename
    std::string new_directory = "../results";
    if (!std::filesystem::exists(new_directory)) {
        std::filesystem::create_directory(new_directory);
    }
    
    std::string base_filename = std::filesystem::path(motion_profile_filename_in).filename().string();
    std::string datetime = getCurrentDateTime();
    std::string filename_without_extension = base_filename.substr(0, base_filename.find_last_of('.'));
    std::string extension = std::filesystem::path(base_filename).extension().string();
    std::string new_base_filename = filename_without_extension + "_ins_ecef" /*+ "_" + datetime*/ + extension;
    std::string motion_profile_filename_out = new_directory + "/" + new_base_filename;

    // Errors filename
    std::string errors_filename_out = new_directory + "/" + 
                                    filename_without_extension + "_ins_ecef_errors" /*+ "_" + datetime*/ + extension;

    // Errors filename
    std::string errors_sigmas_ecef_filename_out = new_directory + "/" + 
                                    filename_without_extension + "_lc_ecef_errors_sigma_ecef" /*+ "_" + datetime*/ + extension;


    // ============== Random gen ==============

    // random device class instance, source of randomness for initializing random seed
    std::random_device rd;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd());

    // ============== IMU parameters ==============

    // Default constructor: tactical grade IMU
    ImuErrors imu_errors;

    // ============== KF config ==============

    KfConfig lc_kf_config;
    // Initial attitude uncertainty per axis (deg, converted to rad)
    lc_kf_config.init_att_unc = deg_to_rad * 1.0;
    // Initial velocity uncertainty per axis (m/s)
    lc_kf_config.init_vel_unc = 0.1;
    // Initial position uncertainty per axis (m)
    lc_kf_config.init_pos_unc = 10.0;
    // Initial accelerometer bias uncertainty per instrument (micro-g, converted
    // to m/s^2)
    lc_kf_config.init_b_a_unc = 1000.0 * micro_g_to_meters_per_second_squared;
    // Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
    lc_kf_config.init_b_g_unc = 10.0 * deg_to_rad / 3600.0;

    // Gyro noise PSD (deg^2 per hour, converted to rad^2/s)                
    lc_kf_config.gyro_noise_PSD = pow(0.02 * deg_to_rad / 60.0, 2.0);
    // Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
    lc_kf_config.accel_noise_PSD = pow(200.0 * micro_g_to_meters_per_second_squared, 2.0);
    // Accelerometer bias random walk PSD (m^2 s^-5)
    lc_kf_config.accel_bias_PSD = 1.0E-7;
    // Gyro bias random walk PSD (rad^2 s^-3)
    lc_kf_config.gyro_bias_PSD = 2.0E-12;

    // Ground truth nav solution in ned
    NavSolutionNed true_nav_ned;
    NavSolutionNed true_nav_ned_old;

    // Ground truth nav solution in ecef
    NavSolutionEcef true_nav_ecef;
    NavSolutionEcef true_nav_ecef_old;

    // Estimated nav solution in ecef
    NavSolutionEcef est_nav_ecef;
    NavSolutionEcef est_nav_ecef_old;

    // Estimated nav solution in ned
    NavSolutionNed est_nav_ned;

    // Ground truth imu measurements from kinematics
    ImuMeasurements true_imu_meas;
    // Simulated imu measurements
    ImuMeasurements imu_meas;
    // Old simulated imu measurements
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();

    // Current time - last time
    // In real use, need all sensors to be synced
    double tor_i;

    // Init motion profile reader & writer
    MotionProfileReader reader(motion_profile_filename_in);
    MotionProfileWriter writer(motion_profile_filename_out);
    ErrorsWriter errors_writer(errors_filename_out);
    ErrorsSigmasEcefWriter errors_sigmas_ecef_writer(errors_sigmas_ecef_filename_out);

    // Error state uncertainty
    Eigen::Matrix<double,15,15> P_matrix = initializeLcPMmatrix(lc_kf_config);

    // Init both true old nav sol, and estimated old nav sol
    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);
    est_nav_ecef_old = true_nav_ecef_old; // Should instead use gnss for this

    while (reader.readNextRow(true_nav_ned)) {

        tor_i = true_nav_ned.time - true_nav_ned_old.time;
        
        // Get true ecef nav sol from motion profile in ned
        true_nav_ecef = nedToEcef(true_nav_ned);

        // ========== IMU SIMULATION ==========

        // Get true specific force and angular rates
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);
        
        // Get imu measurements by applying IMU model
        imu_meas = imuModel(true_imu_meas, imu_meas_old, imu_errors, tor_i, gen);
        // imu_meas = true_imu_meas;

        // ========== NAV EQUATIONS ==========

        // Predict ecef nav solution (INS)

        est_nav_ecef = navEquationsEcef(est_nav_ecef_old, imu_meas, tor_i);
        est_nav_ned = ecefToNed(est_nav_ecef);

        // Compute errors
        ErrorsNed errors = calculateErrorsNed(true_nav_ned, est_nav_ned);

        ErrorsSigmasEcef errors_sigmas_ecef; 
        errors_sigmas_ecef.time = true_nav_ned.time;
        errors_sigmas_ecef.delta_r_eb_e = est_nav_ecef.r_eb_e - true_nav_ecef.r_eb_e;
        errors_sigmas_ecef.delta_v_eb_e = est_nav_ecef.v_eb_e - true_nav_ecef.v_eb_e;
        // errors_sigmas_ecef.delta_eul_eb_e = rToRpy(true_nav_ecef.C_b_e.transpose() * est_nav_ecef.C_b_e);
        errors_sigmas_ecef.delta_eul_eb_e = deSkew(est_nav_ecef.C_b_e * true_nav_ecef.C_b_e.transpose() - Eigen::Matrix3d::Identity());
        
        errors_sigmas_ecef.sigma_delta_r_eb_e << sqrt(P_matrix(6,6)), sqrt(P_matrix(7,7)), sqrt(P_matrix(8,8));
        errors_sigmas_ecef.sigma_delta_v_eb_e << sqrt(P_matrix(3,3)), sqrt(P_matrix(4,4)), sqrt(P_matrix(5,5));
        errors_sigmas_ecef.sigma_delta_eul_eb_e << sqrt(P_matrix(0,0)), sqrt(P_matrix(1,1)), sqrt(P_matrix(2,2));

        // ========== PROPAGATE UNCERTAINTY ==========

        P_matrix  = lcPropUnc(P_matrix, 
                            est_nav_ecef,
                            est_nav_ned,
                            imu_meas,
                            lc_kf_config,
                            tor_i);

        // ========== SAVE RESULTS ==========

        // Save ned output profile
        writer.writeNextRow(est_nav_ned);

        // save errors
        errors_writer.writeNextRow(errors);

        // save errors + sigmas ecef
        errors_sigmas_ecef_writer.writeNextRow(errors_sigmas_ecef);

        true_nav_ecef_old = true_nav_ecef;
        est_nav_ecef_old = est_nav_ecef;
        true_nav_ned_old = true_nav_ned;
        imu_meas_old = imu_meas;

    }

    return 0;
}