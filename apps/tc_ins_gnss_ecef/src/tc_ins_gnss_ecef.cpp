#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>

#include "intnavlib/intnavlib.h"

// To ensure maximum precision, at the cost of speed
// #define EIGEN_DONT_VECTORIZE
// #define EIGEN_FAST_MATH 0

using namespace intnavlib;

int main(int argc, char** argv)
{   

    if(argc != 2) {
        throw std::runtime_error("Pass the motion profile path");
    }

    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

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
    std::string new_base_filename = filename_without_extension + "_tc_ins_gnss_ecef" /*+ "_" + datetime*/ + extension;
    std::string motion_profile_filename_out = new_directory + "/" + new_base_filename;

    // Errors filename
    std::string errors_filename_out = new_directory + "/" + 
                                    filename_without_extension + "_tc_ins_gnss_ecef_errors" /*+ "_" + datetime*/ + extension;
    // Errors filename
    std::string errors_sigmas_ecef_filename_out = new_directory + "/" + 
                                    filename_without_extension + "_tc_ins_gnss_ecef_errors_sigma_ecef" /*+ "_" + datetime*/ + extension;

    // Init motion profile reader & writer
    MotionProfileReader reader(motion_profile_filename_in);
    MotionProfileWriter writer(motion_profile_filename_out);
    ErrorsWriter errors_writer(errors_filename_out);
    ErrorsSigmasEcefWriter errors_sigmas_ecef_writer(errors_sigmas_ecef_filename_out);

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

    GnssConfig gnss_config;

    // Interval between GNSS epochs (s)
    gnss_config.epoch_interval = 0.5;

    // Initial estimated position (m; ECEF)
    gnss_config.init_est_r_ea_e = Eigen::Vector3d::Zero();

    // Number of satellites in constellation
    gnss_config.no_sat = 30.0;
    // Orbital radius of satellites (m)
    gnss_config.r_os = 2.656175E7;
    // Inclination angle of satellites (deg)
    gnss_config.inclination = 55.0;
    // Longitude offset of constellation (deg)
    gnss_config.const_delta_lambda = 0.0;
    // Timing offset of constellation (s)
    gnss_config.const_delta_t = 0.0;

    // Mask angle (deg)
    gnss_config.mask_angle = 0.0; //10.0;
    // Signal in space error SD (m) *Give residual where corrections are applied
    gnss_config.SIS_err_SD = 1.0;
    // Zenith ionosphere error SD (m) *Give residual where corrections are applied
    gnss_config.zenith_iono_err_SD = 2.0;
    // Zenith troposphere error SD (m) *Give residual where corrections are applied
    gnss_config.zenith_trop_err_SD = 0.2;
    // Code tracking error SD (m) *Can extend to account for multipath
    gnss_config.code_track_err_SD = 1.0;
    // Range rate tracking error SD (m/s) *Can extend to account for multipath
    gnss_config.rate_track_err_SD = 0.02;
    // Receiver clock offset at time=0 (m);
    gnss_config.rx_clock_offset = 10000.0;
    // Receiver clock drift at time=0 (m/s);
    gnss_config.rx_clock_drift = 100.0;

    // SD for lc integration
    gnss_config.lc_pos_sd = 2.5;
    gnss_config.lc_vel_sd = 0.1;
    // SDs for tc integration
    gnss_config.pseudo_range_sd = 2.5;
    gnss_config.range_rate_sd = 0.1;

    // ============== KF config ==============

    KfConfig tc_kf_config;
    // Initial attitude uncertainty per axis (deg, converted to rad)
    tc_kf_config.init_att_unc = deg_to_rad * 1.0;
    // Initial velocity uncertainty per axis (m/s)
    tc_kf_config.init_vel_unc = 0.1;
    // Initial position uncertainty per axis (m)
    tc_kf_config.init_pos_unc = 10.0;
    // Initial accelerometer bias uncertainty per instrument (micro-g, converted
    // to m/s^2)
    tc_kf_config.init_b_a_unc = 1000.0 * micro_g_to_meters_per_second_squared;
    // Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
    tc_kf_config.init_b_g_unc = 10.0 * deg_to_rad / 3600.0;
    // Initial clock offset uncertainty per axis (m)
    tc_kf_config.init_clock_offset_unc = 10;
    // Initial clock drift uncertainty per axis (m/s)
    tc_kf_config.init_clock_drift_unc = 0.1;

    // Gyro noise PSD (deg^2 per hour, converted to rad^2/s)                
    tc_kf_config.gyro_noise_PSD = pow(0.02 * deg_to_rad / 60.0, 2.0);
    // Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
    tc_kf_config.accel_noise_PSD = pow(200.0 * micro_g_to_meters_per_second_squared, 2.0);
    // Accelerometer bias random walk PSD (m^2 s^-5)
    tc_kf_config.accel_bias_PSD = 1.0E-7;
    // Gyro bias random walk PSD (rad^2 s^-3)
    tc_kf_config.gyro_bias_PSD = 2.0E-12;
    // Receiver clock frequency-drift PSD (m^2/s^3)
    tc_kf_config.clock_freq_PSD =1;
    // Receiver clock phase-drift PSD (m^2/s)
    tc_kf_config.clock_phase_PSD = 1;

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
    // Estimated clock offset + drift
    double est_clock_offset = 0;
    double est_clock_drift = 0;

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
    Eigen::Matrix<double,17,17> P_matrix = InitializeTcPMmatrix(tc_kf_config);

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
    
    // Estimate initial range biases
    GnssMeasurements gnss_meas_0 = generateGNSSMeasurements(true_nav_ned_old.time,
                                                                sat_pos_vel_0,
                                                                true_nav_ned_old,
                                                                true_nav_ecef_old,
                                                                gnss_biases, 
                                                                gnss_config,
                                                                gen);

    GnssLsPosVelClock gnss_pos_vel_clock_est_0 = gnssLsPositionVelocityClock(gnss_meas_0,
                                                                        est_nav_ecef.r_eb_e,
                                                                        est_nav_ecef.v_eb_e);
    est_clock_offset = gnss_pos_vel_clock_est_0.clock(0);
    est_clock_drift = gnss_pos_vel_clock_est_0.clock(1);

    // Times
    // Current time - last time
    double tor_i;
    double time_last_gnss = 0.0;

    while (reader.readNextRow(true_nav_ned)) {

        auto start = std::chrono::high_resolution_clock::now();

        tor_i = true_nav_ned.time - true_nav_ned_old.time;
        
        // Get true ecef nav sol from motion profile in ned
        true_nav_ecef = nedToEcef(true_nav_ned);

        // ========== IMU SIMULATION ==========

        auto start_imu_sim = std::chrono::high_resolution_clock::now();

        // Get true specific force and angular rates
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);
        // Get imu measurements by applying IMU model
        imu_meas = imuModel(true_imu_meas, imu_meas_old, imu_errors, tor_i, gen);
        // correct imu bias using previous state estimation
        imu_meas.f -= est_acc_bias;
        imu_meas.omega -= est_gyro_bias;

        auto end_imu_sim = std::chrono::high_resolution_clock::now();

        LOG(INFO) << "Imu sim: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_imu_sim - start_imu_sim).count() << "ns";


        // ========== NAV EQUATIONS ==========

        auto start_nav_eqns = std::chrono::high_resolution_clock::now();

        // Predict ecef nav solution (INS)
        est_nav_ecef = navEquationsEcef(est_nav_ecef, imu_meas, tor_i);
        est_nav_ned = ecefToNed(est_nav_ecef);

        auto end_nav_eqns = std::chrono::high_resolution_clock::now();

        LOG(INFO) << "Nav equations: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_nav_eqns - start_nav_eqns).count() << "ns";


        // ========== PROP UNCERTAINTIES ==========

        auto start_unc_prop = std::chrono::high_resolution_clock::now();
            
        // Groves actually puts this in the update part, with a large dt. 
        // But the the underlying approximation hold only
        // if dt is low enough. Therefore slower but better to put it in the prop stage.
        P_matrix  = tcPropUnc(P_matrix, 
                            est_nav_ecef,
                            est_nav_ned,
                            imu_meas,
                            tc_kf_config,
                            tor_i);

        auto end_unc_prop = std::chrono::high_resolution_clock::now();

        LOG(INFO) << "P propagation: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_unc_prop - start_unc_prop).count() << "ns";

        // ========== INTEGRATE POS MEASUREMENTS ==========

        
        double tor_s = true_nav_ned.time - time_last_gnss;
        if( tor_s >= gnss_config.epoch_interval) {

            time_last_gnss = true_nav_ned.time;

            auto start_gnss_sim = std::chrono::high_resolution_clock::now();

            // Simulate GNSS measurement
            SatPosVel sat_pos_vel = satellitePositionsAndVelocities(true_nav_ned.time, gnss_config);
            GnssMeasurements gnss_meas = generateGNSSMeasurements(true_nav_ned.time,
                                                                sat_pos_vel,
                                                                true_nav_ned,
                                                                true_nav_ecef,
                                                                gnss_biases, 
                                                                gnss_config,
                                                                gen);

            auto end_gnss_sim = std::chrono::high_resolution_clock::now();

            LOG(INFO) << "GNSS sim: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_gnss_sim - start_gnss_sim).count() << "ns";

            // KF update -> update posterior
            // if no update, best est is prior
            StateEstEcefTc est_state_ecef_prior;

            est_state_ecef_prior.P_matrix = P_matrix;
            est_state_ecef_prior.nav_sol = est_nav_ecef;
            est_state_ecef_prior.acc_bias = est_acc_bias;
            est_state_ecef_prior.gyro_bias = est_gyro_bias;
            est_state_ecef_prior.clock_offset = est_clock_offset;
            est_state_ecef_prior.clock_drift = est_clock_drift;

            auto start_kf_update = std::chrono::high_resolution_clock::now();

            StateEstEcefTc est_state_ecef_post = tcUpdateKFGnssEcef(gnss_meas, 
                                                                    est_state_ecef_prior,
                                                                    tor_s);

            auto end_kf_update = std::chrono::high_resolution_clock::now();

            LOG(INFO) << "GNSS KF Update: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end_kf_update - start_kf_update).count() << "ns";

            P_matrix = est_state_ecef_post.P_matrix;
            est_nav_ecef = est_state_ecef_post.nav_sol;
            est_acc_bias = est_state_ecef_post.acc_bias;
            est_gyro_bias = est_state_ecef_post.gyro_bias;
            est_clock_offset = est_state_ecef_post.clock_offset;
            est_clock_drift = est_state_ecef_post.clock_drift;

        }

        LOG(INFO) << " ============== ";
        
        // ========== COMPUTE ERRORS ==========

        est_nav_ned = ecefToNed(est_nav_ecef);
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

        // ========== SAVE RESULTS ==========

        // Save ned output profile
        writer.writeNextRow(est_nav_ned);

        // save errors
        errors_writer.writeNextRow(errors);

        // save errors + sigmas ecef
        errors_sigmas_ecef_writer.writeNextRow(errors_sigmas_ecef);

        true_nav_ecef_old = true_nav_ecef;
        true_nav_ned_old = true_nav_ned;
        imu_meas_old = imu_meas;
    }

    return 0;
}