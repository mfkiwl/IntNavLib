#ifndef HELPERS_H
#define HELPERS_H

#include <Eigen/Dense>

namespace helpers {

    // Constants for WGS84 model
    constexpr double R_0 = 6378137.0;  // WGS84 Equatorial radius in meters
    constexpr double e = 0.0818191908425; // WGS84 eccentricity
    constexpr double omega_ie = 7.292115e-5;  // Earth rotation rate
    constexpr double mu = 3.986004418e14; // WGS84 Earth gravitational constant (m^3 s^-2)
    constexpr double J_2 = 1.082627e-3; // WGS84 Earth's second gravitational constant
    
    // IMU measurements
    struct ImuMeasurements {
        Eigen::Vector3d f; // Specific forces
        Eigen::Vector3d omega; // angular velocities
    };

    // IMU errors
    struct ImuErrors {

    };

    struct NavSolutionNed {
        double time;
        double latitude;  // Latitude in radians
        double longitude; // Longitude in radians
        double height;    // Height in meters
        Eigen::Vector3d v_b_n; // Velocity in NED frame (North, East, Down)
        Eigen::Matrix3d C_b_n; // Body-to-NED rotation matrix
    };

    struct NavSolutionEcef {
        double time;
        Eigen::Vector3d p_b_e; // Position in ECEF frame
        Eigen::Vector3d v_b_e; // Velocity in ECEF frame
        Eigen::Matrix3d C_b_e; // Body-to-ECEF rotation matrix
    };

    // IMU quantization residuals
    typedef std::array<double, 6> ImuQuantResiduals;

    // Nav equations in ECEF
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav_sol, 
                                const ImuMeasurements & imu_meas);

    // Gets real imu measurements from kinematics
    ImuMeasurements kinematicsEcef(const NavSolutionEcef & old_nav, 
                                    const NavSolutionEcef & new_nav);
    
    // Gets gravity vector at given ecef position, in ecef frame
    Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e);

    // IMU model
    ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                            const ImuErrors & imu_errors);

    // NED to ECEF
    NavSolutionEcef nedToEcef(const NavSolutionNed & nav_sol_ned);

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


};

#endif