#ifndef SIMULATION_H
#define SIMULATION_H

#include <eigen3/Eigen/Dense>

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    // Gets real imu measurements from kinematics
    ImuMeasurements kinematicsEcef(const NavSolutionEcef & new_nav,
                                    const NavSolutionEcef & old_nav);

    // IMU model
    ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                            const ImuMeasurements & old_imu_meas,
                            const ImuErrors & imu_errors,
                            const double & tor_i,
                            std::mt19937 & gen);

    // Generic ecef position sensor model
    PosMeasEcef genericPosSensModel(const NavSolutionEcef & true_nav, 
                                const double & pos_sigma,
                                std::mt19937 & gen);

    // Generic ecef pose sensor model
    PosRotMeasEcef genericPosRotSensModel(const NavSolutionEcef & true_nav, 
                                        const double & pos_sigma,
                                        const double & rot_sigma,
                                        std::mt19937 & gen);

    // Generate satellites positions and velocities given config. returns ECEF Cartesian positions and
    // ECEF velocities of all satellites in the constellation. Simple circular
    // orbits with regularly distributed satellites are modeled.
    SatPosVel satellitePositionsAndVelocities(const double & time, const GnssConfig & GNSS_config);

    // GNSS receiver sensor model: Generates a set of pseudo-range and pseudo-
    // range rate measurements for all satellites above the elevation mask angle 
    // and adds satellite positions and velocities to the datesets.
    GnssMeasurements generateGNSSMeasurements(const double & time,
                                        const SatPosVel & gnss_pos_vel,
                                        const NavSolutionNed& true_nav_ned,
                                        const NavSolutionEcef& true_nav_ecef,
                                        const Eigen::VectorXd& GNSS_biases, 
                                        const GnssConfig& GNSS_config,
                                        std::mt19937 & gen);

    // Initializes the GNSS range errors due to signal
    // in space, ionosphere and troposphere errors based on the elevation angles. 
    Eigen::VectorXd initializeGNSSBiases(const NavSolutionEcef & true_nav_ecef,
                                    const NavSolutionNed & true_nav_ned,
                                    const SatPosVel & gnss_pos_vel,
                                    const GnssConfig& GNSS_config,
                                    std::mt19937 & gen);

};

#endif