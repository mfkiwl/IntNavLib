#include <thread>
#include <chrono>
#include <random>

#include "helpers.h"
#include "lc_ins_gnss.h"
#include "motion_profile_reader.h"

int main(int argc, char** argv)
{   

    if(argc != 2) {
        throw std::runtime_error("Pass the motion profile path");
    }

    std::string motion_profile_filename(argv[1]);

    MotionProfileReader reader(motion_profile_filename);

    // ============== Random gen ==============

    // random device class instance, source of 'true' randomness for initializing random seed
    std::random_device rd; 

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd()); 

    // ============== IMU parameters ==============

    ImuErrors imu_errors;

    // Accelerometer biases (micro-g, converted to m/s^2; body axes)
    imu_errors.b_a << 900,-1300,800;
    imu_errors.b_a = imu_errors.b_a * micro_g_to_meters_per_second_squared;

    // Gyro biases (deg/hour, converted to rad/sec; body axes)
    imu_errors.b_g << -9, 13, -8;
    imu_errors.b_g = imu_errors.b_g * deg_to_rad / 3600;

    // Accelerometer scale factor and cross coupling errors (ppm, converted to
    // unitless; body axes)
    imu_errors.M_a << 500, -300, 200,
                    -150, -600, 250,
                    -250,  100, 450;
    imu_errors.M_a = imu_errors.M_a * 1e-6;

    // Gyro scale factor and cross coupling errors (ppm, converted to unitless;
    // body axes)
    imu_errors.M_g << 400, -300,  250,
                        0, -300, -150,
                        0,    0, -350; 
    imu_errors.M_g = imu_errors.M_g * 1e-6;

    // Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes)
    imu_errors.G_g << 0.9, -1.1, -0.6,
                    -0.5,  1.9, -1.6,
                    0.3,  1.1, -1.3;
    imu_errors.G_g = imu_errors.G_g * deg_to_rad / (3600 * 9.80665);  

    // Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5)                
    imu_errors.accel_noise_root_PSD = 100 * micro_g_to_meters_per_second_squared;

    // Gyro noise root PSD (deg per root hour, converted to rad s^-0.5)                
    imu_errors.gyro_noise_root_PSD = 0.01 * deg_to_rad / 60;

    // Accelerometer quantization level (m/s^2)
    imu_errors.accel_quant_level = 1e-2;

    // Gyro quantization level (rad/s)
    imu_errors.gyro_quant_level = 2e-4;

    // ============== 

    // Ground truth nav solution in ned
    NavSolutionNed true_nav_ned;
    NavSolutionNed true_nav_ned_old;

    // Ground truth nav solution in ecef
    NavSolutionEcef true_nav_ecef;
    NavSolutionEcef true_nav_ecef_old;

    // Estimated nav solution in ecef
    NavSolutionEcef est_nav_ecef;
    NavSolutionEcef est_nav_ecef_old;

    // Ground truth imu measurements from kinematics
    ImuMeasurements true_imu_meas;

    // Simulated imu measurements
    ImuMeasurements imu_meas;

    // Current time - last time
    // In real use, need all sensors to be synced
    double tor_i;

    // Init both true old nav sol, and estimated old nav sol
    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);
    est_nav_ecef_old = true_nav_ecef_old; // Should instead use gnss for this

    while (reader.readNextRow(true_nav_ned)) {

        tor_i = true_nav_ned.time - true_nav_ned_old.time;
        
        // Get true ecef nav sol from motion profile in ned
        true_nav_ecef = nedToEcef(true_nav_ned);

        // TODO unit test

        // Get true specific force and angular rates
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // TODO unit test
        
        // Get imu measurements by applying IMU model
        imu_meas = imuModel(true_imu_meas, imu_errors, tor_i, gen);

        // Predict ecef nav solution (INS)
        est_nav_ecef = navEquationsEcef(est_nav_ecef_old, imu_meas, tor_i);

        // TODO Plot true nav sol & est nav sol on map

        true_nav_ecef_old = true_nav_ecef;
        est_nav_ecef_old = est_nav_ecef;

        true_nav_ned_old = true_nav_ned;

    }

    return 0;
}