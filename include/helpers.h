#ifndef HELPERS_H
#define HELPERS_H

#include <Eigen/Dense>

namespace helpers {

    // Constants for WGS84 model
    constexpr double R_0 = 6378137.0;  // WGS84 Equatorial radius in meters
    constexpr double e = 0.0818191908425; // WGS84 eccentricity

    // IMU measurements
    struct ImuMeasurements {
        
    };

    // IMU errors
    struct ImuErrors {

    };

    struct NavSolutionNed {
        double time;
        double latitude;  // Latitude in radians
        double longitude; // Longitude in radians
        double height;    // Height in meters
        Eigen::Vector3d velocity_ned; // Velocity in NED frame (North, East, Down)
        Eigen::Matrix3d C_b_n; // Body-to-NED rotation matrix
    };

    struct NavSolutionEcef {
        double time;
        Eigen::Vector3d position_ecef; // Position in ECEF frame
        Eigen::Vector3d velocity_ecef; // Velocity in ECEF frame
        Eigen::Matrix3d C_b_e;         // Body-to-ECEF rotation matrix
    };

    // IMU quantization residuals
    typedef std::array<double, 6> ImuQuantResiduals;

    // Nav equations in ECEF
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav_sol, 
                                const ImuMeasurements & imu_meas);

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

    // Gets real imu measurements from kinematics
    ImuMeasurements kinematicsEcef(const NavSolutionEcef & old_nav, const NavSolutionEcef & new_nav);

    Eigen::Matrix3d rpyToR(const double & roll, const double & pitch, const double & yaw);

};

#endif