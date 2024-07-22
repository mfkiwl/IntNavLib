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
    
};

#endif