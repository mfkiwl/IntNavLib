#ifndef HELPERS_H
#define HELPERS_H

#include <random>
#include <ctime>
#include <iomanip>
#include <iostream>

#include <eigen3/Eigen/Dense>

#include "constants_types.h"

namespace intnavlib {

    /// \defgroup Helpers
    /// @{
    
    /// \brief Gets gravity vector at given ECEF position, in ECEF frame.
    /// \param[in] r_eb_e Position vector in ECEF frame.
    /// \return Gravity vector in ECEF frame.
    Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e);

    /// \brief Converts a navigation solution from NED (North, East, Down) frame to ECEF (Earth-Centered, Earth-Fixed) frame.
    /// \param[in] nav_sol_ned Navigation solution in NED frame.
    /// \return Navigation solution in ECEF frame.
    NavSolutionEcef nedToEcef(const NavSolutionNed & nav_sol_ned);

    /// \brief Converts a navigation solution from ECEF (Earth-Centered, Earth-Fixed) frame to NED (North, East, Down) frame.
    /// \param[in] nav_sol_ecef Navigation solution in ECEF frame.
    /// \return Navigation solution in NED frame.
    NavSolutionNed ecefToNed(const NavSolutionEcef & nav_sol_ecef);

    /// \brief Initializes the error covariance matrix for a navigation Kalman Filter.
    /// \param[in] tc_kf_config Configuration parameters for the Kalman Filter,
    /// including initial uncertainties for attitude, velocity, position, biases,
    /// and clock offset/drift.
    /// \return An Eigen::Matrix<double,17,17> representing the initialized error covariance matrix.
    Eigen::Matrix<double,17,17> initializePMmatrix(const KfConfig & tc_kf_config);

    /// \brief Init navigation filter state from ground truth.
    /// \param[in] true_nav_ecef Ground truth navigation solution.
    /// \param[in] kf_config Navigation filter configuration.
    /// \param[in] gnss_meas GNSS measurements.
    /// \param[in] gen Random source.
    /// \return Navigation filter state.
    StateEstEcef initStateFromGroundTruth(const NavSolutionEcef & true_nav_ecef, const KfConfig & kf_config, const GnssMeasurements & gnss_meas, std::mt19937 & gen);

    /// \brief Compute errors from state estimate and ground truth.
    /// \param[in] state_est_ecef State estimate.
    /// \param[in] true_nav_ecef grounud truth.
    /// \return Errors and standard deviations.
    ErrorsSigmasEcef getErrorsSigmasEcef(const StateEstEcef & state_est_ecef, const NavSolutionEcef & true_nav_ecef);

    /// \brief Converts degrees to radians.
    /// \param[in] degrees Angle in degrees.
    /// \return Angle in radians.
    inline double degToRad(const double & degrees) {
        return degrees * 0.01745329252;
    }

    /// \brief Converts radians to degrees.
    /// \param[in] rads Angle in radians.
    /// \return Angle in degrees.
    inline double radToDeg(const double & rads) {
        return rads / 0.01745329252;
    }

    /// \brief Converts Euler angles (roll, pitch, yaw) to a rotation matrix.
    /// \param[in] rpy Eigen::Vector3d containing roll, pitch, and yaw angles (in radians).
    /// \return 3x3 rotation matrix.
    Eigen::Matrix3d rpyToR(const Eigen::Vector3d & rpy);

    /// \brief Converts a rotation matrix to Euler angles (roll, pitch, yaw).
    /// \param[in] C 3x3 rotation matrix.
    /// \return Eigen::Vector3d containing roll, pitch, and yaw angles (in radians).
    Eigen::Vector3d rToRpy(const Eigen::Matrix3d & C);

    /// \brief Gets the skew-symmetric matrix from a 3D vector.
    /// For a vector `a = [ax, ay, az]`, the skew-symmetric matrix `S` is:
    /// `[ 0  -az  ay ]`
    /// `[ az  0  -ax ]`
    /// `[-ay  ax  0  ]`
    /// \param[in] a Input 3D vector.
    /// \return 3x3 skew-symmetric matrix.
    Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d & a);

    /// \brief Signum function.
    /// Returns -1 if `val < 0`, 0 if `val == 0`, and 1 if `val > 0`.
    /// \tparam T Type of the input value.
    /// \param[in] val Input value.
    /// \return Sign of the input value.
    template <typename T> inline T sgn(const T & val) {
    return T((T(0) < val) - (val < T(0)));
    }

    /// \brief Gets the current date and time as a formatted string.
    /// The format is "YYYY-MM-DD_HH-MM-SS".
    /// \return A string representing the current date and time.
    std::string getCurrentDateTime();

    /// \brief Gets the radii of curvature (meridian and transverse) from latitude.
    /// \param[in] L Latitude in radians.
    /// \return Eigen::Vector2d where x is the meridian radius and y is the transverse radius.
    Eigen::Vector2d radiiOfCurvature(double L);

    /// \brief Calculates position, velocity, and attitude errors in the local NED (North, East, Down) frame.
    /// \param[in] true_nav_sol True navigation solution in NED frame.
    /// \param[in] est_nav_sol Estimated navigation solution in NED frame.
    /// \return Structure containing the calculated errors.
    ErrorsNed calculateErrorsNed(const NavSolutionNed & true_nav_sol, 
                                const NavSolutionNed & est_nav_sol);

    /// \brief Extracts the vector from a skew-symmetric matrix.
    /// This is the inverse operation of `skewSymmetric`.
    /// For a skew-symmetric matrix `S`:
    /// `[ 0  -az  ay ]`
    /// `[ az  0  -ax ]`
    /// `[-ay  ax  0  ]`
    /// the function returns the vector `[ax, ay, az]`.
    /// \param[in] S 3x3 skew-symmetric matrix.
    /// \return The 3D vector from which the skew-symmetric matrix was formed.
    Eigen::Vector3d deSkew(const Eigen::Matrix3d & S);

    /// \brief Calculates position, velocity, clock offset, and clock drift rate using unweighted iterated least squares.
    /// \param[in] GNSS_measurements Structure containing GNSS pseudo-range and pseudo-range rate measurements,
    ///                              along with satellite positions and velocities.
    /// \param[in] prior_r_ea_e Prior estimated antenna position in ECEF frame.
    /// \param[in] prior_v_ea_e Prior estimated antenna velocity in ECEF frame.
    /// \return Structure containing the estimated antenna position, velocity, clock offset and drift rate.
    GnssLsPosVelClock gnssLsPositionVelocityClock(const GnssMeasurements & GNSS_measurements,
                                                const Eigen::Vector3d & prior_r_ea_e,
                                                const Eigen::Vector3d & prior_v_ea_e);
    
    /// \brief Calculates position, velocity using unweighted iterated least squares.
    /// \param[in] GNSS_measurements Structure containing GNSS pseudo-range and pseudo-range rate measurements,
    ///                              along with satellite positions and velocities.
    /// \param[in] prior_r_ea_e Prior estimated antenna position in ECEF frame.
    /// \param[in] prior_v_ea_e Prior estimated antenna velocity in ECEF frame.
    /// \param[in] gnss_config GNSS config parameters
    /// \return Structure containing the estimated antenna position, velocity.
    GnssPosVelMeasEcef gnssLsPositionVelocity(const GnssMeasurements & GNSS_measurements,
                                    const Eigen::Vector3d & prior_r_ea_e,
                                    const Eigen::Vector3d & prior_v_ea_e,
                                    const GnssConfig & gnss_config);
    
    /// \brief Get IMU errors for a tactical grade IMU.
    ImuErrors tacticalImuErrors();
    /// \brief Get default GNSS config.
    GnssConfig defaultGnssConfig();
    /// \brief Get default KfConfig.
    KfConfig tacticalImuKFConfig();                  
                                
    /// @}                      
                                
};

#endif