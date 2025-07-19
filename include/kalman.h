#ifndef KALMAN_H
#define KALMAN_H

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

/// \brief Initializes the error covariance matrix for a Loosely Coupled (LC) Kalman Filter.
/// \param[in] lc_kf_config Configuration parameters for the LC Kalman Filter,
/// including initial uncertainties for attitude, velocity, position, and biases.
/// \return An Eigen::Matrix<double,15,15> representing the initialized error covariance matrix.
Eigen::Matrix<double,15,15> initializeLcPMmatrix(const KfConfig & lc_kf_config);

/// \brief Initializes the error covariance matrix for a Tightly Coupled (TC) Kalman Filter.
/// \param[in] tc_kf_config Configuration parameters for the TC Kalman Filter,
/// including initial uncertainties for attitude, velocity, position, biases,
/// and clock offset/drift.
/// \return An Eigen::Matrix<double,17,17> representing the initialized error covariance matrix.
Eigen::Matrix<double,17,17> initializeTcPMmatrix(const KfConfig & tc_kf_config);

/// \brief Propagates the error-state uncertainty for a Loosely Coupled (LC) Kalman Filter.
/// This function updates the error covariance matrix based on the system dynamics
/// and noise characteristics over a given integration time.
/// \param[in] P_matrix_old The error covariance matrix from the previous time step.
/// \param[in] old_nav_est_ecef The estimated navigation solution in ECEF frame at the previous time step.
/// \param[in] old_nav_est_ned The estimated navigation solution in NED frame at the previous time step.
/// \param[in] imu_meas IMU measurements used for propagation.
/// \param[in] lc_kf_config Configuration parameters for the LC Kalman Filter,
/// including noise PSDs for gyro, accelerometer, and biases.
/// \param[in] tor_s Integration time (delta_t) in seconds.
/// \return An Eigen::Matrix<double,15,15> representing the propagated error covariance matrix.
Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & lc_kf_config,
                                        const double & tor_s);
/// \brief Propagates the error-state uncertainty for a Tightly Coupled (TC) Kalman Filter.
/// This function updates the error covariance matrix based on the system dynamics
/// and noise characteristics over a given integration time, including clock states.
/// \param[in] P_matrix_old The error covariance matrix from the previous time step.
/// \param[in] old_nav_est_ecef The estimated navigation solution in ECEF frame at the previous time step.
/// \param[in] old_nav_est_ned The estimated navigation solution in NED frame at the previous time step.
/// \param[in] imu_meas IMU measurements used for propagation.
/// \param[in] tc_kf_config Configuration parameters for the TC Kalman Filter,
/// including noise PSDs for gyro, accelerometer, biases, and clock states.
/// \param[in] tor_s Integration time (delta_t) in seconds.
/// \return An Eigen::Matrix<double,17,17> representing the propagated error covariance matrix.
Eigen::Matrix<double,17,17> tcPropUnc(const Eigen::Matrix<double,17,17> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & tc_kf_config,
                                        const double & tor_s);

/// \brief Performs a Loosely Coupled (LC) Kalman Filter update using an ECEF position measurement.
/// This is a closed-loop error-state KF update that corrects the navigation solution
/// and IMU bias estimates by integrating a position measurement.
/// \param[in] pos_meas The position measurement in ECEF frame.
/// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
/// and its covariance matrix.
/// \return A StateEstEcefLc structure containing the updated state estimation.
StateEstEcefLc lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                    const StateEstEcefLc & state_est_old);

/// \brief Performs a Loosely Coupled (LC) Kalman Filter update using a GNSS position and velocity measurement.
/// This is a closed-loop error-state KF update that corrects the navigation solution
/// and IMU bias estimates by integrating a GNSS position and velocity measurement.
/// \param[in] pos_vel_gnss_meas The GNSS position and velocity measurement in ECEF frame.
/// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
/// and its covariance matrix.
/// \return A StateEstEcefLc structure containing the updated state estimation.
StateEstEcefLc lcUpdateKFGnssEcef (const GnssPosVelMeasEcef & pos_vel_gnss_meas, 
                                    const StateEstEcefLc & state_est_old);

/// \brief Performs a Tightly Coupled (TC) Kalman Filter update using GNSS pseudo-range and pseudo-range rate measurements.
/// This is a closed-loop error-state KF update that corrects the navigation solution,
/// IMU bias estimates, and receiver clock states by integrating raw GNSS measurements.
/// \param[in] gnss_meas The GNSS measurements, including pseudo-ranges, pseudo-range rates,
/// satellite positions, and satellite velocities.
/// \param[in] state_est_prior The prior estimated state (navigation solution, biases, and clock states)
/// and its covariance matrix.
/// \param[in] tor_s Elapsed time since last update in seconds.
/// \return A StateEstEcefTc structure containing the updated state estimation.
StateEstEcefTc tcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcefTc & state_est_prior,
                                    const double & tor_s);

/// \brief Performs a Loosely Coupled (LC) Kalman Filter update using an ECEF position and rotation measurement.
/// This is a closed-loop error-state KF update that corrects the navigation solution
/// and IMU bias estimates by integrating a position and rotation measurement.
/// \param[in] pos_rot_meas The position and rotation measurement in ECEF frame.
/// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
/// and its covariance matrix.
/// \return A StateEstEcefLc structure containing the updated state estimation.
StateEstEcefLc lcUpdateKFPosRotEcef (const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcefLc & state_est_old);

};



#endif