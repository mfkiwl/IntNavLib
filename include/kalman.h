#ifndef KALMAN_H
#define KALMAN_H

#include "constants_types.h"

namespace intnavlib {

Eigen::Matrix<double,15,15> InitializeLcPMmatrix(const LcKfConfig & lc_kf_config) {

Eigen::Matrix<double,15,15> p_matrix;

// Initialize error covariance matrix
p_matrix = Eigen::Matrix<double,15,15>::Zero();
p_matrix.block<3,3>(0,0) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_att_unc,2); // attitude error
p_matrix.block<3,3>(3,3) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_vel_unc,2); // vel error
p_matrix.block<3,3>(6,6) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_pos_unc,2); // pos error
p_matrix.block<3,3>(9,9) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_b_a_unc,2); // acc bias error
p_matrix.block<3,3>(12,12) = Eigen::Matrix3d::Identity() * pow(lc_kf_config.init_b_g_unc,2); // gyro bias error

return p_matrix;
}

};

#endif