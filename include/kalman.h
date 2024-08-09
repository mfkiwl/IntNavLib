#ifndef KALMAN_H
#define KALMAN_H

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

// Initialize error covariance matrix
Eigen::Matrix<double,15,15> InitializeLcPMmatrix(const KfConfig & lc_kf_config);

// Initialize error covariance matrix
Eigen::Matrix<double,17,17> InitializeTcPMmatrix(const KfConfig & tc_kf_config);

// Propagates error-state uncertainty, given old P and integration time
Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & lc_kf_config,
                                        const double & tor_s);

// Propagates error-state uncertainty, given old P and integration time
Eigen::Matrix<double,17,17> tcPropUnc(const Eigen::Matrix<double,17,17> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & lc_kf_config,
                                        const double & tor_s);

// Loosely coupled, ECEF, closed loop error-state KF update
// Updates error state and thus navigation solution + bias estimate, integrating a position measurement
StateEstEcefLc lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                    const StateEstEcefLc & state_est_old);

// Loosely coupled, ECEF, closed loop error-state KF update
// Updates error state and thus navigation solution + bias estimate, integrating a gnss position + vel measurement
StateEstEcefLc lcUpdateKFGnssEcef (const GnssPosVelMeasEcef & pos_vel_gnss_meas, 
                                    const StateEstEcefLc & state_est_old);

// Tightly coupled, ECEF, closed loop error-state KF update
// Updates error state and thus navigation solution + bias estimate, integrating gnss ranges and range rates + sat pos & vel
StateEstEcefTc tcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcefTc & state_est_prior,
                                    const double & tor_s);

// Loosely coupled, ECEF, closed loop error-state KF update
// Updates error state and thus navigation solution + bias estimate, integrating a position + rotation measurement
StateEstEcefLc lcUpdateKFPosRotEcef (const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcefLc & state_est_old);

};



#endif