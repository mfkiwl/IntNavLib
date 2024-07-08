#ifndef HELPERS_H
#define HELPERS_H

#include <random>
#include <ctime>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>


namespace helpers {

    // Epsilon for == 0 checks
    // Todo tune
    constexpr double EPSILON = 1.0e-100; 

    // Constants for WGS84 model
    constexpr double R_0 = 6378137.0;  // WGS84 Equatorial radius in meters
    constexpr double e = 0.0818191908425; // WGS84 eccentricity
    constexpr double omega_ie = 7.292115e-5;  // Earth rotation rate
    constexpr double mu = 3.986004418e14; // WGS84 Earth gravitational constant (m^3 s^-2)
    constexpr double J_2 = 1.082627e-3; // WGS84 Earth's second gravitational constant

    constexpr double deg_to_rad = 0.01745329252;
    constexpr double rad_to_deg = 1.0/deg_to_rad;
    constexpr double micro_g_to_meters_per_second_squared = 9.80665e-6;
    
    // IMU measurements
    struct ImuMeasurements {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        double time;
        Eigen::Vector3d f; // Specific forces
        Eigen::Vector3d omega; // angular velocities
        Eigen::Vector3d quant_residuals_f;
        Eigen::Vector3d quant_residuals_omega;
    };

    // IMU errors
    struct ImuErrors {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Eigen::Vector3d b_a;              // Accelerometer biases (m/s^2)
        Eigen::Vector3d b_g;              // Gyro biases (rad/s)
        Eigen::Matrix3d M_a;              // Accelerometer scale factor and cross coupling errors
        Eigen::Matrix3d M_g;              // Gyro scale factor and cross coupling errors            
        Eigen::Matrix3d G_g;              // Gyro g-dependent biases (rad-sec/m)             
        double accel_noise_root_PSD;   // Accelerometer noise root PSD (m s^-1.5)
        double gyro_noise_root_PSD;    // Gyro noise root PSD (rad s^-0.5)
        double accel_quant_level;      // Accelerometer quantization level (m/s^2)
        double gyro_quant_level;       // Gyro quantization level (rad/s)
    };

    struct NavSolutionNed {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        double time;
        double latitude;  // Latitude in radians
        double longitude; // Longitude in radians
        double height;    // Height in meters
        Eigen::Vector3d v_eb_n; // Velocity in NED frame (North, East, Down)
        Eigen::Matrix3d C_b_n; // Body-to-NED rotation matrix
    };

    struct NavSolutionEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        double time;
        Eigen::Vector3d r_eb_e; // Position in ECEF frame
        Eigen::Vector3d v_eb_e; // Velocity in ECEF frame
        Eigen::Matrix3d C_b_e; // Body-to-ECEF rotation matrix
    };

    struct ErrorsNed {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        double time;
        Eigen::Vector3d delta_r_eb_n; // NED Position error
        Eigen::Vector3d delta_v_eb_n; // NED Velocity error
        Eigen::Vector3d delta_eul_nb_n; // NED orientation error
    };

    // Nav equations in ECEF
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i);

    // Gets real imu measurements from kinematics
    ImuMeasurements kinematicsEcef(const NavSolutionEcef & new_nav,
                                    const NavSolutionEcef & old_nav);
    
    // Gets gravity vector at given ecef position, in ecef frame
    Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e);

    // IMU model
    ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                            const ImuErrors & imu_errors,
                            const double & tor_i,
                            std::mt19937 & gen);

    // NED to ECEF
    NavSolutionEcef nedToEcef(const NavSolutionNed & nav_sol_ned);

    // ECEF to NED
    NavSolutionNed ecefToNed(const NavSolutionEcef & nav_sol_ecef);

    // Converts degrees to radians
    inline double degToRad(const double & degrees) {
        return degrees * 0.01745329252;
    }

    // Converts radians to degrees
    inline double radToDeg(const double & rads) {
        return rads / 0.01745329252;
    }

    Eigen::Matrix3d rpyToR(const Eigen::Vector3d & rpy);

    Eigen::Vector3d rToRpy(const Eigen::Matrix3d & C);

    Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d & a);

    // signum function
    template <typename T> inline T sgn(const T & val) {
    return T((T(0) < val) - (val < T(0)));
    }

    // Function to get the current date and time as a formatted string
    std::string getCurrentDateTime();

    // Get radii of curvature from latitude
    Eigen::Vector2d radiiOfCurvature(double L);

    ErrorsNed calculateErrorsNed(const NavSolutionNed & true_nav_sol, 
                                const NavSolutionNed & est_nav_sol);

};

#endif