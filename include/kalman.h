#ifndef KALMAN_H
#define KALMAN_H

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

// Initialize error covariance matrix
Eigen::Matrix<double,15,15> InitializeLcPMmatrix(const LcKfConfig & lc_kf_config);

// Propagates error-state uncertainty, given old P and integration time
Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const LcKfConfig & lc_kf_config,
                                        const double & tor_s);

};

#endif