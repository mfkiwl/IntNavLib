#include <chrono>
#include <memory>
#include <queue>
#include <thread>
#include <mutex>

#include "rclcpp/rclcpp.hpp"
#include "sensor_msgs/msg/imu.hpp"
#include "sensor_msgs/msg/nav_sat_fix.hpp"
#include "geometry_msgs/msg/pose_with_covariance_stamped.hpp"
#include "intnavlib/intnavlib.h"

#define MAX_BUFFER_SIZE_IMU 200
#define MAX_BUFFER_SIZE_GNSS 50
#define TIME_EPSILON_S 0.01

using namespace std::chrono_literals;
using namespace intnavlib;

class InsGnssNode : public rclcpp::Node {

public:

    InsGnssNode() : Node("ins_gnss_node") {

        // Initialize variables from ROS2 parameters
        init();

        // Initialize subscribers and publisher
        imu_sub_ = this->create_subscription<sensor_msgs::msg::Imu>(
            imu_topic_, 100, std::bind(&InsGnssNode::imuCallback, this, std::placeholders::_1));
        gnss_sub_ = this->create_subscription<sensor_msgs::msg::NavSatFix>(
            gnss_topic_, 10, std::bind(&InsGnssNode::gnssCallback, this, std::placeholders::_1));
        pose_pub_ = this->create_publisher<geometry_msgs::msg::PoseWithCovarianceStamped>(out_pose_topic_, 100);

        // Start processing thread
        processing_thread_ = std::thread(&InsGnssNode::processData, this);
    }

    ~InsGnssNode() {
        processing_thread_.join();
    }

private:

    std::string imu_topic_;
    std::string gnss_topic_;
    std::string out_pose_topic_;

    rclcpp::Subscription<sensor_msgs::msg::Imu>::SharedPtr imu_sub_;
    rclcpp::Subscription<sensor_msgs::msg::NavSatFix>::SharedPtr gnss_sub_;
    rclcpp::Publisher<geometry_msgs::msg::PoseWithCovarianceStamped>::SharedPtr pose_pub_;

    // Buffers
    std::queue<sensor_msgs::msg::Imu::SharedPtr> imu_buffer_;
    std::queue<sensor_msgs::msg::NavSatFix::SharedPtr> gnss_buffer_;
    std::mutex imu_buffer_mutex_;
    std::condition_variable imu_cv;
    std::mutex gnss_buffer_mutex_;

    // Init params
    std::vector<double> init_lla_;
    std::vector<double> init_rpy_b_n_;
    std::vector<double> init_v_eb_n_;

    // Navigation state
    NavSolutionEcef est_nav_ecef_;
    Eigen::Vector3d est_acc_bias_;
    Eigen::Vector3d est_gyro_bias_;
    Eigen::Matrix<double, 15, 15> P_matrix_;

    // EKF config
    KfConfig lc_kf_config_;

    std::thread processing_thread_;

    // If log, write results to file
    bool write_log_;
    std::unique_ptr<MotionProfileWriter> log_writer_;

    void init() {

        // Declare parameters

        // Topics
        declare_parameter<std::string>("imu_topic");
        declare_parameter<std::string>("gnss_topic");
        declare_parameter<std::string>("out_pose_topic");

        // Logging
        declare_parameter<bool>("write_log");
        declare_parameter<std::string>("log_path");

        // EKF config
        declare_parameter<double>("kf_config.init_att_unc");
        declare_parameter<double>("kf_config.init_vel_unc");
        declare_parameter<double>("kf_config.init_pos_unc");
        declare_parameter<double>("kf_config.init_b_a_unc");
        declare_parameter<double>("kf_config.init_b_g_unc");
        declare_parameter<double>("kf_config.gyro_noise_PSD");
        declare_parameter<double>("kf_config.accel_noise_PSD");
        declare_parameter<double>("kf_config.accel_bias_PSD");
        declare_parameter<double>("kf_config.gyro_bias_PSD");

        // Init 
        declare_parameter<std::vector<double>>("init_lla");
        declare_parameter<std::vector<double>>("init_v_eb_n");
        declare_parameter<std::vector<double>>("init_rpy_b_n");

        // Get parameters

        // Topics 
        imu_topic_ = get_parameter("imu_topic").as_string();
        gnss_topic_ = get_parameter("gnss_topic").as_string();
        out_pose_topic_ = get_parameter("out_pose_topic").as_string();

        // Logging
        write_log_ = get_parameter("write_log").as_bool();
        std::string log_path = get_parameter("log_path").as_string();
        if(write_log_){
            log_writer_ = std::make_unique<MotionProfileWriter>(log_path);
        }

        // EKF config
        lc_kf_config_.init_att_unc = get_parameter("kf_config.init_att_unc").as_double();
        lc_kf_config_.init_vel_unc = get_parameter("kf_config.init_vel_unc").as_double();
        lc_kf_config_.init_pos_unc = get_parameter("kf_config.init_pos_unc").as_double();
        lc_kf_config_.init_b_a_unc = get_parameter("kf_config.init_b_a_unc").as_double();
        lc_kf_config_.init_b_g_unc = get_parameter("kf_config.init_b_g_unc").as_double();
        lc_kf_config_.gyro_noise_PSD = get_parameter("kf_config.gyro_noise_PSD").as_double();
        lc_kf_config_.accel_noise_PSD = get_parameter("kf_config.accel_noise_PSD").as_double();
        lc_kf_config_.accel_bias_PSD = get_parameter("kf_config.accel_bias_PSD").as_double();
        lc_kf_config_.gyro_bias_PSD = get_parameter("kf_config.gyro_bias_PSD").as_double();

        // Init nav solution
        init_lla_ = get_parameter("init_lla").as_double_array();
        init_v_eb_n_ = get_parameter("init_v_eb_n").as_double_array();
        init_rpy_b_n_ = get_parameter("init_rpy_b_n").as_double_array();

        NavSolutionNed est_nav_ned = NavSolutionNed{0.0,
                                                    deg_to_rad * init_lla_[0], 
                                                    deg_to_rad * init_lla_[1], 
                                                    init_lla_[2], 
                                                    Eigen::Vector3d(init_v_eb_n_[0], init_v_eb_n_[1], init_v_eb_n_[2]), 
                                                    rpyToR(deg_to_rad * Eigen::Vector3d(init_rpy_b_n_[0], init_rpy_b_n_[1], init_rpy_b_n_[2])).transpose()};
        est_nav_ecef_ = nedToEcef(est_nav_ned);

        est_acc_bias_ = Eigen::Vector3d::Zero();
        est_gyro_bias_ = Eigen::Vector3d::Zero();
        P_matrix_ = InitializeLcPMmatrix(lc_kf_config_);

        // Log loaded parameters for debugging
        
        RCLCPP_INFO(this->get_logger(), "IMU topic: %s", imu_topic_.c_str());
        RCLCPP_INFO(this->get_logger(), "GNSS topic: %s", gnss_topic_.c_str());

        RCLCPP_INFO(this->get_logger(), "Init_lla = %f %f %f", init_lla_[0], init_lla_[1], init_lla_[2]);
        RCLCPP_INFO(this->get_logger(), "Init_v_eb_n = %f %f %f", init_v_eb_n_[0], init_v_eb_n_[1], init_v_eb_n_[2]);
        RCLCPP_INFO(this->get_logger(), "Init_rpy_b_n = %f %f %f", init_rpy_b_n_[0], init_rpy_b_n_[1], init_rpy_b_n_[2]);

        RCLCPP_INFO(this->get_logger(), "Init_att_unc = %f", lc_kf_config_.init_att_unc);
        RCLCPP_INFO(this->get_logger(), "Init_vel_unc = %f", lc_kf_config_.init_vel_unc);
        RCLCPP_INFO(this->get_logger(), "Init_pos_unc = %f", lc_kf_config_.init_pos_unc);

        RCLCPP_INFO(this->get_logger(), "Nav Filter Initialized! Waiting on measurements...");
    }

    void imuCallback(const sensor_msgs::msg::Imu::SharedPtr msg) {
        //  Threadsafe push to buffer, notify condition variable
        std::lock_guard<std::mutex> lock(imu_buffer_mutex_);
        if(imu_buffer_.size() < MAX_BUFFER_SIZE_IMU) {
            imu_buffer_.push(msg);
            imu_cv.notify_one();
        }
        else RCLCPP_WARN(this->get_logger(), "IMU buffer full: dropping!");
    }

    void gnssCallback(const sensor_msgs::msg::NavSatFix::SharedPtr msg) {
        //  Threadsafe push to buffer
        std::lock_guard<std::mutex> lock(gnss_buffer_mutex_);
        if(gnss_buffer_.size() < MAX_BUFFER_SIZE_GNSS)
            gnss_buffer_.push(msg);
        else RCLCPP_WARN(this->get_logger(), "GNSS buffer full: dropping!");
    }

    void processData() {
        while (rclcpp::ok()) {

            sensor_msgs::msg::Imu::SharedPtr imu_msg;
            sensor_msgs::msg::NavSatFix::SharedPtr gnss_msg;

            // Threadsafe read from message buffers
            // Want to get oldest IMU measurement (front)
            double imu_time = 0;
            {   
                // Wait on IMU buffer
                std::unique_lock<std::mutex> lock(imu_buffer_mutex_);
                imu_cv.wait(lock, [this] { return !imu_buffer_.empty();});

                // Get measurement timestamp
                rclcpp::Time imu_stamp(imu_buffer_.front()->header.stamp);
                imu_time = imu_stamp.seconds();
                
                // If new, use measurement, else just pop
                if (imu_time > est_nav_ecef_.time)
                    imu_msg = imu_buffer_.front();
                imu_buffer_.pop();
                
            }

            // If msg received but not good (old), continue
            if(!imu_msg) continue;

            // We received a new IMU measurement: do propagation

            RCLCPP_INFO(this->get_logger(), "Predict: %f", imu_time);

            ImuMeasurements imu_meas;
            imu_meas.f = Eigen::Vector3d(
                imu_msg->linear_acceleration.x,
                imu_msg->linear_acceleration.y,
                imu_msg->linear_acceleration.z);
            imu_meas.omega = Eigen::Vector3d(
                imu_msg->angular_velocity.x,
                imu_msg->angular_velocity.y,
                imu_msg->angular_velocity.z);
            imu_meas.time = imu_time;

            double tor_i = imu_meas.time - est_nav_ecef_.time;

            // Apply bias corrections
            imu_meas.f -= est_acc_bias_;
            imu_meas.omega -= est_gyro_bias_;

            // Predict
            est_nav_ecef_ = navEquationsEcef(est_nav_ecef_, imu_meas, tor_i);

            // Propagate uncertainties
            // Consider doing this at GNSS rate, if need to speed up
            {
                NavSolutionNed est_nav_ned = ecefToNed(est_nav_ecef_);
                P_matrix_ = lcPropUnc(P_matrix_, est_nav_ecef_, est_nav_ned, imu_meas, lc_kf_config_, tor_i);
            }

            // Now, see if we also have a new GNSS measurement

            {
                std::lock_guard<std::mutex> lock(gnss_buffer_mutex_);

                if (!gnss_buffer_.empty()) {
                    rclcpp::Time gnss_stamp(gnss_buffer_.front()->header.stamp);
                    // RCLCPP_INFO(this->get_logger(), "GNSS time %f", gnss_stamp.seconds());
                    // RCLCPP_INFO(this->get_logger(), "IMU time %f", imu_time);
                    if (abs(gnss_stamp.seconds() - imu_time) < TIME_EPSILON_S)
                        gnss_msg = gnss_buffer_.front();
                    gnss_buffer_.pop();
                    
                }
            }

            // If new GNSS, do update

            if (gnss_msg) {

                rclcpp::Time gnss_stamp(gnss_msg->header.stamp);
                double gnss_time = gnss_stamp.seconds();
                RCLCPP_INFO(this->get_logger(), "GNSS Update: %f", gnss_time);

                PosMeasEcef gnss_meas;
                gnss_meas.time = gnss_time;
                gnss_meas.r_eb_e = nedToEcef(NavSolutionNed{0,gnss_msg->latitude, gnss_msg->longitude, gnss_msg->altitude, Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity()}).r_eb_e;
                gnss_meas.cov_mat = Eigen::Matrix3d::Identity() * std::pow(2.5,2); // Covariance matrix

                // Update navigation state using GNSS
                StateEstEcefLc est_state_ecef_prior{est_nav_ecef_, est_acc_bias_, est_gyro_bias_, P_matrix_};
                StateEstEcefLc est_state_ecef_post = lcUpdateKFPosEcef(gnss_meas, est_state_ecef_prior);

                // Update state variables
                est_nav_ecef_ = est_state_ecef_post.nav_sol;
                est_acc_bias_ = est_state_ecef_post.acc_bias;
                est_gyro_bias_ = est_state_ecef_post.gyro_bias;
                P_matrix_ = est_state_ecef_post.P_matrix;
            }

            // Publish estimated pose
            publishPose();

            // Log nav sol to file
            if(write_log_){
                NavSolutionNed est_nav_ned = ecefToNed(est_nav_ecef_);
                log_writer_->writeNextRow(est_nav_ned);
            }
        }
    }

    void publishPose() {
        // Create a PoseWithCovarianceStamped message
        auto pose_msg = geometry_msgs::msg::PoseWithCovarianceStamped();
        pose_msg.header.stamp = rclcpp::Time(static_cast<int64_t>(est_nav_ecef_.time * 1e9));
        pose_msg.header.frame_id = "ecef";

        // Fill in position
        pose_msg.pose.pose.position.x = est_nav_ecef_.r_eb_e[0];
        pose_msg.pose.pose.position.y = est_nav_ecef_.r_eb_e[1];
        pose_msg.pose.pose.position.z = est_nav_ecef_.r_eb_e[2];

        // Fill in orientation (convert rotation matrix to quaternion)
        Eigen::Quaterniond q(est_nav_ecef_.C_b_e);
        pose_msg.pose.pose.orientation.x = q.x();
        pose_msg.pose.pose.orientation.y = q.y();
        pose_msg.pose.pose.orientation.z = q.z();
        pose_msg.pose.pose.orientation.w = q.w();

        Eigen::Matrix<double, 6, 6> pose_cov;
        pose_cov.block<3,3>(0,0) = P_matrix_.block<3,3>(6,6);  // position covariance
        pose_cov.block<3,3>(0,3) = P_matrix_.block<3,3>(6,0);  // cross-covariance
        pose_cov.block<3,3>(3,0) = P_matrix_.block<3,3>(0,6);  // cross-covariance
        pose_cov.block<3,3>(3,3) = P_matrix_.block<3,3>(0,0);  // attitude covariance

        // Copy the 6x6 matrix into the covariance array (which is in row-major order)
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                pose_msg.pose.covariance[i * 6 + j] = pose_cov(i, j);
            }
        }

        pose_pub_->publish(pose_msg);
    }
};

int main(int argc, char *argv[]) {
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<InsGnssNode>());
    rclcpp::shutdown();
    return 0;
}
