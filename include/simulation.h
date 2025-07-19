#ifndef SIMULATION_H
#define SIMULATION_H

#include <eigen3/Eigen/Dense>

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    /// \brief Calculates true IMU measurements (specific force and angular velocity) from kinematic states.
    /// This function determines the IMU outputs that would be observed given the change in navigation states.
    /// \param[in] new_nav The navigation solution at the current time step.
    /// \param[in] old_nav The navigation solution at the previous time step.
    /// \return An ImuMeasurements structure containing the true specific force and angular velocity.
    ImuMeasurements kinematicsEcef(const NavSolutionEcef & new_nav,
                                    const NavSolutionEcef & old_nav);

    /// \brief Models IMU sensor outputs by adding errors (biases, noise, quantization) to true IMU measurements.
    /// \param[in] true_imu_meas The true IMU measurements (specific force and angular velocity).
    /// \param[in] old_imu_meas The IMU measurements from the previous time step, used for quantization residuals.
    /// \param[in] imu_errors Structure containing IMU error parameters (biases, scale factors, noise PSDs, quantization levels).
    /// \param[in] tor_i The integration time interval (delta_t) in seconds.
    /// \param[in] gen A Mersenne Twister random number generator for noise generation.
    /// \return An ImuMeasurements structure containing the modeled IMU outputs with errors.
    ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                            const ImuMeasurements & old_imu_meas,
                            const ImuErrors & imu_errors,
                            const double & tor_i,
                            std::mt19937 & gen);

    /// \brief Models a generic ECEF position sensor by adding Gaussian noise to the true position.
    /// \param[in] true_nav The true navigation solution in ECEF frame.
    /// \param[in] pos_sigma The standard deviation of the position measurement noise (in meters).
    /// \param[in] gen A Mersenne Twister random number generator for noise generation.
    /// \return A PosMeasEcef structure containing the noisy position measurement and its covariance.
    PosMeasEcef genericPosSensModel(const NavSolutionEcef & true_nav, 
                                const double & pos_sigma,
                                std::mt19937 & gen);

    /// \brief Models a generic ECEF pose (position and rotation) sensor by adding Gaussian noise to the true position and rotation.
    /// \param[in] true_nav The true navigation solution in ECEF frame.
    /// \param[in] pos_sigma The standard deviation of the position measurement noise (in meters).
    /// \param[in] rot_sigma The standard deviation of the rotation measurement noise (in radians).
    /// \param[in] gen A Mersenne Twister random number generator for noise generation.
    /// \return A PosRotMeasEcef structure containing the noisy position and rotation measurement and its covariance.
    PosRotMeasEcef genericPosRotSensModel(const NavSolutionEcef & true_nav, 
                                        const double & pos_sigma,
                                        const double & rot_sigma,
                                        std::mt19937 & gen);

    /// \brief Generates satellite positions and velocities in the ECEF frame based on a simplified circular orbit model.
    /// \param[in] time Current simulation time in seconds.
    /// \param[in] GNSS_config Configuration parameters for the GNSS constellation.
    /// \return A SatPosVel structure containing the ECEF Cartesian positions and velocities of all satellites.
    SatPosVel satellitePositionsAndVelocities(const double & time, const GnssConfig & GNSS_config);

    /// \brief Generates GNSS pseudo-range and pseudo-range rate measurements for satellites above the elevation mask angle.
    /// This function also includes satellite positions and velocities in the output.
    /// \param[in] time Current simulation time in seconds.
    /// \param[in] gnss_pos_vel Satellite positions and velocities.
    /// \param[in] true_nav_ned True navigation solution in NED frame.
    /// \param[in] true_nav_ecef True navigation solution in ECEF frame.
    /// \param[in] GNSS_biases Pre-calculated GNSS range biases for each satellite.
    /// \param[in] GNSS_config Configuration parameters for the GNSS receiver and constellation.
    /// \param[in] gen A Mersenne Twister random number generator for noise generation.
    /// \return A GnssMeasurements structure containing the generated pseudo-range and pseudo-range rate measurements,
    ///         along with corresponding satellite positions and velocities.
    GnssMeasurements generateGNSSMeasurements(const double & time,
                                        const SatPosVel & gnss_pos_vel,
                                        const NavSolutionNed& true_nav_ned,
                                        const NavSolutionEcef& true_nav_ecef,
                                        const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MAX_GNSS_SATELLITES> & GNSS_biases, 
                                        const GnssConfig& GNSS_config,
                                        std::mt19937 & gen);

    /// \brief Initializes GNSS range biases due to signal-in-space, ionosphere, and troposphere errors.
    /// These biases are calculated based on satellite elevation angles and configured error standard deviations.
    /// \param[in] true_nav_ecef True navigation solution in ECEF frame.
    /// \param[in] true_nav_ned True navigation solution in NED frame.
    /// \param[in] gnss_pos_vel Satellite positions and velocities.
    /// \param[in] GNSS_config Configuration parameters for the GNSS constellation and error models.
    /// \param[in] gen A Mersenne Twister random number generator for noise generation.
    /// \return An Eigen::Matrix containing the initialized GNSS range biases for each satellite.
    Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MAX_GNSS_SATELLITES>
    initializeGNSSBiases(const NavSolutionEcef & true_nav_ecef,
                                    const NavSolutionNed & true_nav_ned,
                                    const SatPosVel & gnss_pos_vel,
                                    const GnssConfig& GNSS_config,
                                    std::mt19937 & gen);

};

#endif