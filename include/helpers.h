#ifndef HELPERS_H
#define HELPERS_H

#include <Eigen/Dense>

namespace helpers {

    // IMU measurements
    struct ImuMeasurements {
        
    };

    // IMU errors
    struct ImuErrors {

    };

    // Navigation solution in ecef
    struct NavSolutionEcef {

    };

    // Navigation solution in ned
    struct NavSolutionNed {

    };

    // Struct to hold a row of the motion profile
    struct MotionProfileRow {
        double time;
        double latitude;
        double longitude;
        double height;
        double north_velocity;
        double east_velocity;
        double down_velocity;
        double roll_angle;
        double pitch_angle;
        double yaw_angle;
    };

    // IMU quantization residuals
    typedef std::array<double, 6> ImuQuantResiduals;

    // Nav equations in ECEF
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav_sol, 
                                const ImuMeasurements & imu_meas);

    // IMU model
    ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                            const ImuErrors & imu_errors
                            );

    // NED to ECEF
    NavSolutionEcef nedToEcef(const NavSolutionNed & nav_sol_ned);

    inline double degToRad(const double & degrees) {
        return degrees * 0.01745329252; // Conversion factor from degrees to radians
    }

}

#endif