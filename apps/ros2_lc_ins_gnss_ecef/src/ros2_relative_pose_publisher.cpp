#include "rclcpp/rclcpp.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "geometry_msgs/msg/transform_stamped.hpp"
#include <eigen3/Eigen/Dense>

class RelativePoseNode : public rclcpp::Node
{
public:
    RelativePoseNode() : Node("relative_pose_node")
    {
        estimated_pose_sub_ = this->create_subscription<geometry_msgs::msg::PoseStamped>(
            "/nav_sol/ecef", 10, std::bind(&RelativePoseNode::estimatedPoseCallback, this, std::placeholders::_1));

        ground_truth_pose_sub_ = this->create_subscription<geometry_msgs::msg::PoseStamped>(
            "/ground_truth/pose", 10, std::bind(&RelativePoseNode::groundTruthPoseCallback, this, std::placeholders::_1));

        relative_pose_pub_ = this->create_publisher<geometry_msgs::msg::PoseStamped>("/relative_pose", 10);
    }

private:
    void estimatedPoseCallback(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
    {
        estimated_pose_ = *msg;
        computeAndPublishRelativePose();
    }

    void groundTruthPoseCallback(const geometry_msgs::msg::PoseStamped::SharedPtr msg)
    {
        ground_truth_pose_ = *msg;
        computeAndPublishRelativePose();
    }

    void computeAndPublishRelativePose()
    {
        if (!estimated_pose_ || !ground_truth_pose_)
            return;

        // Convert ground truth and estimated poses to Eigen objects
        Eigen::Vector3d r_gt(
            ground_truth_pose_->pose.position.x,
            ground_truth_pose_->pose.position.y,
            ground_truth_pose_->pose.position.z);

        Eigen::Vector3d r_est(
            estimated_pose_->pose.position.x,
            estimated_pose_->pose.position.y,
            estimated_pose_->pose.position.z);

        Eigen::Quaterniond q_gt(
            ground_truth_pose_->pose.orientation.w,
            ground_truth_pose_->pose.orientation.x,
            ground_truth_pose_->pose.orientation.y,
            ground_truth_pose_->pose.orientation.z);

        Eigen::Quaterniond q_est(
            estimated_pose_->pose.orientation.w,
            estimated_pose_->pose.orientation.x,
            estimated_pose_->pose.orientation.y,
            estimated_pose_->pose.orientation.z);

        // Compute relative position
        Eigen::Vector3d relative_position = q_gt.inverse() * (r_est - r_gt);

        // Compute relative orientation
        Eigen::Quaterniond relative_orientation = q_gt.inverse() * q_est;

        // Publish relative pose
        auto relative_pose_msg = geometry_msgs::msg::PoseStamped();
        relative_pose_msg.header.stamp = this->now();
        relative_pose_msg.header.frame_id = "world";

        relative_pose_msg.pose.position.x = relative_position.x();
        relative_pose_msg.pose.position.y = relative_position.y();
        relative_pose_msg.pose.position.z = relative_position.z();

        relative_pose_msg.pose.orientation.x = relative_orientation.x();
        relative_pose_msg.pose.orientation.y = relative_orientation.y();
        relative_pose_msg.pose.orientation.z = relative_orientation.z();
        relative_pose_msg.pose.orientation.w = relative_orientation.w();

        relative_pose_pub_->publish(relative_pose_msg);

        RCLCPP_INFO(this->get_logger(), "Published relative pose");
    }

    rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr estimated_pose_sub_;
    rclcpp::Subscription<geometry_msgs::msg::PoseStamped>::SharedPtr ground_truth_pose_sub_;
    rclcpp::Publisher<geometry_msgs::msg::PoseStamped>::SharedPtr relative_pose_pub_;

    std::optional<geometry_msgs::msg::PoseStamped> estimated_pose_;
    std::optional<geometry_msgs::msg::PoseStamped> ground_truth_pose_;
};

int main(int argc, char **argv)
{
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<RelativePoseNode>());
    rclcpp::shutdown();
    return 0;
}
