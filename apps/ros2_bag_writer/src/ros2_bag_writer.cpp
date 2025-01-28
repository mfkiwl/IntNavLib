#include <thread>
#include <chrono>
#include <random>
#include <filesystem>
#include <glog/logging.h>

#include <rclcpp/rclcpp.hpp>
#include <sensor_msgs/msg/imu.hpp>
#include <sensor_msgs/msg/nav_sat_fix.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <rosbag2_cpp/writer.hpp>
#include <rosbag2_cpp/writer_interfaces/base_writer_interface.hpp>

#include "intnavlib/intnavlib.h"

using namespace intnavlib;

int main(int argc, char** argv)
{
    if(argc != 2) {
        throw std::runtime_error("Pass the motion profile path");
    }

    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    // ============== Random gen ==============

    // random device class instance, source of randomness for initializing random seed
    std::random_device rd;

    // Mersenne twister PRNG, initialized with seed from previous random device instance
    std::mt19937 gen(rd());

    // ============== IMU parameters ==============

    ImuErrors imu_errors;

    // Accelerometer biases (micro-g, converted to m/s^2; body axes)gnss_msg
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

    // ============== GNSS parameters ==============

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
    gnss_config.mask_angle = 10.0;
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

    // ==========================================

    // Initialize ROS 2
    rclcpp::init(argc, argv);
    auto node = rclcpp::Node::make_shared("bag_writer_node");

    // Input profile filename
    std::string motion_profile_filename_in(argv[1]);

    // Initialize motion profile reader
    MotionProfileReader reader(motion_profile_filename_in);

    // Initialize ROS 2 bag writer
    auto bag_writer = std::make_shared<rosbag2_cpp::Writer>();
    std::string bag_directory = "output/Profile.bag";
    if (std::filesystem::exists("output/"))
        std::filesystem::remove_all("output/");
    std::filesystem::create_directory("output/");
    bag_writer->open(bag_directory);

    // Declare message objects
    sensor_msgs::msg::Imu imu_clean_msg;
    sensor_msgs::msg::Imu imu_dirty_msg;
    sensor_msgs::msg::NavSatFix gnss_msg;
    geometry_msgs::msg::PoseStamped gt_pose_msg;

    // Persistent variables
    NavSolutionNed true_nav_ned;
    NavSolutionNed true_nav_ned_old;
    NavSolutionEcef true_nav_ecef;
    NavSolutionEcef true_nav_ecef_old;
    ImuMeasurements true_imu_meas;
    ImuMeasurements imu_meas;
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();
    double time_last_gnss = 0.0;

    // Initialize true navigation solution
    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);

    // Init GNSS range biases
    SatPosVel sat_pos_vel_0 = satellitePositionsAndVelocities(true_nav_ned_old.time, 
                                                            gnss_config);
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MAX_GNSS_SATELLITES> gnss_biases = 
                                    initializeGNSSBiases(true_nav_ecef_old,
                                                        true_nav_ned_old,
                                                        sat_pos_vel_0,
                                                        gnss_config,
                                                        gen);

    // gnss_biases *= 0;

    while (reader.readNextRow(true_nav_ned)) {

        LOG(INFO) << "Writing t = " << true_nav_ned.time;

        double tor_i = true_nav_ned.time - true_nav_ned_old.time;

        // Update true navigation solution in ECEF
        true_nav_ecef = nedToEcef(true_nav_ned);

        // Generate true IMU measurements
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // Populate IMU message
        imu_clean_msg.header.stamp = rclcpp::Time(static_cast<int64_t>(true_nav_ned.time * 1e9));
        imu_clean_msg.linear_acceleration.x = true_imu_meas.f[0];
        imu_clean_msg.linear_acceleration.y = true_imu_meas.f[1];
        imu_clean_msg.linear_acceleration.z = true_imu_meas.f[2];
        imu_clean_msg.angular_velocity.x = true_imu_meas.omega[0];
        imu_clean_msg.angular_velocity.y = true_imu_meas.omega[1];
        imu_clean_msg.angular_velocity.z = true_imu_meas.omega[2];
        bag_writer->write(imu_clean_msg, "imu_clean/data", imu_clean_msg.header.stamp);

        // Get imu measurements by applying IMU model
        imu_meas = imuModel(true_imu_meas, imu_meas_old, imu_errors, tor_i, gen);
        // Populate IMU message
        imu_dirty_msg.header.stamp = imu_clean_msg.header.stamp;
        imu_dirty_msg.linear_acceleration.x = imu_meas.f[0];
        imu_dirty_msg.linear_acceleration.y = imu_meas.f[1];
        imu_dirty_msg.linear_acceleration.z = imu_meas.f[2];
        imu_dirty_msg.angular_velocity.x = imu_meas.omega[0];
        imu_dirty_msg.angular_velocity.y = imu_meas.omega[1];
        imu_dirty_msg.angular_velocity.z = imu_meas.omega[2];
        bag_writer->write(imu_dirty_msg, "imu_dirty/data", imu_clean_msg.header.stamp);

        // Populate Pose message
        gt_pose_msg.header.stamp = imu_clean_msg.header.stamp;
        gt_pose_msg.header.frame_id = "world";
        gt_pose_msg.pose.position.x = true_nav_ecef.r_eb_e[0];
        gt_pose_msg.pose.position.y = true_nav_ecef.r_eb_e[1];
        gt_pose_msg.pose.position.z = true_nav_ecef.r_eb_e[2];

        Eigen::Quaterniond q(true_nav_ecef.C_b_e);
        gt_pose_msg.pose.orientation.x = q.x();
        gt_pose_msg.pose.orientation.y = q.y();
        gt_pose_msg.pose.orientation.z = q.z();
        gt_pose_msg.pose.orientation.w = q.w();
        bag_writer->write(gt_pose_msg, "ground_truth/pose", imu_clean_msg.header.stamp);

        // Every n secs, write gnss msg
        double tor_s = true_nav_ned.time - time_last_gnss;
        if( tor_s >= gnss_config.epoch_interval) {

            LOG(INFO) << "Writing GNSS";

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
            
            // Estimate receiver pos + vel with NL LS
            GnssLsPosVelClock pos_vel_clock_gnss_meas_ecef = gnssLsPositionVelocityClock(gnss_meas,
                                                                        true_nav_ecef.r_eb_e,
                                                                        true_nav_ecef.v_eb_e);

            // Convert position to LLA
            NavSolutionNed lla_gnss_meas = ecefToNed(NavSolutionEcef{0,pos_vel_clock_gnss_meas_ecef.r_ea_e, Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity()});

            // Populate GNSS message
            gnss_msg.header.stamp = imu_clean_msg.header.stamp;
            gnss_msg.latitude = lla_gnss_meas.latitude;
            gnss_msg.longitude = lla_gnss_meas.longitude;
            gnss_msg.altitude = lla_gnss_meas.height;
            bag_writer->write(gnss_msg, "gnss/fix", imu_clean_msg.header.stamp);

        }

        true_nav_ecef_old = true_nav_ecef;
        true_nav_ned_old = true_nav_ned;
        imu_meas_old = imu_meas;
    }

    bag_writer->close();
    rclcpp::shutdown();
    return 0;
}
