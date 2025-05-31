from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node

def generate_launch_description():

    log_level_arg = DeclareLaunchArgument(
        'log_level',
        default_value='INFO',
        description='Logging level (e.g., DEBUG, INFO, WARN, ERROR, FATAL)'
    )

    log_level = LaunchConfiguration('log_level')

    main_node = Node(
        package='ros2_lc_ins_gnss_ecef',
        executable='ros2_lc_ins_gnss_ecef',
        name='ros2_lc_ins_gnss_ecef',
        output='screen',
        parameters=['config/Profile_3_config.yaml'],
        arguments=['--ros-args', '--log-level', log_level],
    )

    rviz_node = Node(
        package='rviz2',
        executable='rviz2',
        name='rviz2',
        arguments=['-d', 'config/config.rviz'],
        output='screen'
    )

    static_pub_node = Node(
        package="tf2_ros",
        executable="static_transform_publisher",
        name="static_tf_enu_to_profile_3",
        arguments=["4069880.333114", "-255758.024327", "4900749.559590",
                   "0.246771", "0.231751", "0.644153", "0.685903", "ecef", "enu"],
        output="screen",
    )

    return LaunchDescription([
        log_level_arg,
        main_node,
        rviz_node,
        static_pub_node
    ])
