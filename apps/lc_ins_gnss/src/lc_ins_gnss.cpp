#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include "intnavlib/intnavlib.h"

// To ensure full precisiion
// #define EIGEN_DONT_VECTORIZE
// #define EIGEN_FAST_MATH 0

using namespace intnavlib;

int main(int argc, char** argv)
{   

    if(argc != 2) {
        throw std::runtime_error("Pass the motion profile path");
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
    std::string new_base_filename = filename_without_extension + "_ins" /*+ "_" + datetime*/ + extension;
    std::string motion_profile_filename_out = new_directory + "/" + new_base_filename;

    // Errors filename
    std::string errors_filename_out = new_directory + "/" + 
                                    filename_without_extension + "_ins_errors" /*+ "_" + datetime*/ + extension;

    // ============== Random gen ==============

    // random device class instance, source of randomness for initializing random seed
    std::random_device rd;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd());

    // ============== IMU parameters ==============

    ImuErrors imu_errors;

    // Accelerometer biases (micro-g, converted to m/s^2; body axes)
    imu_errors.b_a << 900.0,-1300.0,800.0;
    imu_errors.b_a = imu_errors.b_a * micro_g_to_meters_per_second_squared;

    // Gyro biases (deg/hour, converted to rad/sec; body axes)
    imu_errors.b_g << -9.0, 13.0, -8.0;
    imu_errors.b_g = imu_errors.b_g * deg_to_rad / 3600.0;

    // Accelerometer scale factor and cross coupling errors (ppm, converted to
    // unitless; body axes)
    imu_errors.M_a << 500.0, -300.0, 200.0,
                    -150.0, -600.0, 250.0,
                    -250.0,  100.0, 450.0;
    imu_errors.M_a = imu_errors.M_a * 1.0e-6;

    // Gyro scale factor and cross coupling errors (ppm, converted to unitless;
    // body axes)
    imu_errors.M_g << 400.0, -300.0,  250.0,
                        0.0, -300.0, -150.0,
                        0.0,    0.0, -350.0; 
    imu_errors.M_g = imu_errors.M_g * 1.0e-6;

    // Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes)
    imu_errors.G_g << 0.9, -1.1, -0.6,
                    -0.5,  1.9, -1.6,
                    0.3,  1.1, -1.3;
    imu_errors.G_g = imu_errors.G_g * deg_to_rad / (3600.0 * 9.80665);  

    // Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5)                
    imu_errors.accel_noise_root_PSD = 100.0 * micro_g_to_meters_per_second_squared;

    // Gyro noise root PSD (deg per root hour, converted to rad s^-0.5)                
    imu_errors.gyro_noise_root_PSD = 0.01 * deg_to_rad / 60.0;

    // Accelerometer quantization level (m/s^2)
    imu_errors.accel_quant_level = 1.0e-2;

    // Gyro quantization level (rad/s)
    imu_errors.gyro_quant_level = 2.0e-4;

    // ============== GNSS config ==============

    GnssConfig gnssConfig;

    // Interval between GNSS epochs (s)
    gnssConfig.epoch_interval = 0.5;

    // Initial estimated position (m; ECEF)
    gnssConfig.init_est_r_ea_e = Eigen::Vector3d::Zero();

    // Number of satellites in constellation
    gnssConfig.no_sat = 30.0;
    // Orbital radius of satellites (m)
    gnssConfig.r_os = 2.656175E7;
    // Inclination angle of satellites (deg)
    gnssConfig.inclination = 55.0;
    // Longitude offset of constellation (deg)
    gnssConfig.const_delta_lambda = 0.0;
    // Timing offset of constellation (s)
    gnssConfig.const_delta_t = 0.0;

    // Mask angle (deg)
    gnssConfig.mask_angle = 10.0;
    // Signal in space error SD (m) *Give residual where corrections are applied
    gnssConfig.SIS_err_SD = 1.0;
    // Zenith ionosphere error SD (m) *Give residual where corrections are applied
    gnssConfig.zenith_iono_err_SD = 2.0;
    // Zenith troposphere error SD (m) *Give residual where corrections are applied
    gnssConfig.zenith_trop_err_SD = 0.2;
    // Code tracking error SD (m) *Can extend to account for multipath
    gnssConfig.code_track_err_SD = 1.0;
    // Range rate tracking error SD (m/s) *Can extend to account for multipath
    gnssConfig.rate_track_err_SD = 0.02;
    // Receiver clock offset at time=0 (m);
    gnssConfig.rx_clock_offset = 10000.0;
    // Receiver clock drift at time=0 (m/s);
    gnssConfig.rx_clock_drift = 100.0;

    // ============== KF config ==============

    LcKfConfig kf_config;

    // Initial attitude uncertainty per axis (deg, converted to rad)
    kf_config.init_att_unc = deg_to_rad * 1.0;
    // Initial velocity uncertainty per axis (m/s)
    kf_config.init_vel_unc = 0.1;
    // Initial position uncertainty per axis (m)
    kf_config.init_pos_unc = 10.0;
    // Initial accelerometer bias uncertainty per instrument (micro-g, converted
    // to m/s^2)
    kf_config.init_b_a_unc = 1000.0 * micro_g_to_meters_per_second_squared;
    // Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
    kf_config.init_b_g_unc = 10.0 * deg_to_rad / 3600.0;

    // Gyro noise PSD (deg^2 per hour, converted to rad^2/s)                
    kf_config.gyro_noise_PSD = pow(0.02 * deg_to_rad / 60.0, 2.0);
    // Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
    kf_config.accel_noise_PSD = pow(200.0 * micro_g_to_meters_per_second_squared, 2.0);
    // Accelerometer bias random walk PSD (m^2 s^-5)
    kf_config.accel_bias_PSD = 1.0E-7;
    // Gyro bias random walk PSD (rad^2 s^-3)
    kf_config.gyro_bias_PSD = 2.0E-12;

    // ============

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

    // Current time - last time
    // In real use, need all sensors to be synced
    double tor_i;

    // Init motion profile reader & writer
    MotionProfileReader reader(motion_profile_filename_in);
    MotionProfileWriter writer(motion_profile_filename_out);
    ErrorsWriter errors_writer(errors_filename_out);

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
        imu_meas = imuModel(true_imu_meas, imu_errors, tor_i, gen);

        // ========== NAV EQUATIONS ==========

        // Predict ecef nav solution (INS)
        est_nav_ecef = navEquationsEcef(est_nav_ecef_old, imu_meas, tor_i);
        est_nav_ned = ecefToNed(est_nav_ecef);

        // Compute errors
        ErrorsNed errors = calculateErrorsNed(true_nav_ned, est_nav_ned);

        // INTEGRATE POS MEASUREMENTS
        // if(true_nav_ned.time - time_last_gnss> )

        // ========== SAVE RESULTS ==========

        // Save ned output profile
        writer.writeNextRow(est_nav_ned);

        // save errors
        errors_writer.writeNextRow(errors);

        true_nav_ecef_old = true_nav_ecef;
        est_nav_ecef_old = est_nav_ecef;
        true_nav_ned_old = true_nav_ned;

    }

    return 0;
}