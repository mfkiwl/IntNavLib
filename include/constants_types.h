
#ifndef CONSTANTS_TYPES_H
#define CONSTANTS_TYPES_H

#include <eigen3/Eigen/Dense>

namespace intnavlib {

/// \defgroup constants_types Constants and Types
/// @{

/// Max number of GNSS satellites.
/// Leaving some extra room for future expansion.
constexpr int kMaxGnssSatellites = 35;
/// Epsilon for checking if a value is approximately zero.
constexpr double kEpsilon = 1.0e-100; 
/// WGS84 Equatorial radius in meters.
constexpr double kR0 = 6378137.0;  // WGS84 Equatorial radius in meters
/// WGS84 eccentricity.
constexpr double kEccentricity = 0.0818191908425; // WGS84 eccentricity
/// Earth rotation rate.
constexpr double kOmega_ie = 7.292115e-5;  // Earth rotation rate
/// WGS84 Earth gravitational constant (m^3 s^-2).
constexpr double kGravConst = 3.986004418e14; // WGS84 Earth gravitational constant (m^3 s^-2)
/// WGS84 Earth's second gravitational constant.
constexpr double kJ2 = 1.082627e-3; // WGS84 Earth's second gravitational constant
/// Speed of light in m/s.
constexpr double c = 299792458.0 ; // Speed of light in m/s
/// Conversion factor from degrees to radians.
constexpr double kDegToRad = 0.01745329252;
/// Conversion factor from radians to degrees.
constexpr double kRadToDeg = 1.0/kDegToRad;
/// Conversion factor from micro-g to meters per second squared.
constexpr double kMuGToMetersPerSecondSquared = 9.80665e-6;

/// \brief Structure to hold IMU measurements.
struct ImuMeasurements {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the measurement.
    double time;
    /// Specific forces (accelerometer measurements).
    Eigen::Vector3d f;
    /// Angular velocities (gyroscope measurements).
    Eigen::Vector3d omega;
    /// Quantization residuals for specific forces.
    Eigen::Vector3d quant_residuals_f;
    /// Quantization residuals for angular velocities.
    Eigen::Vector3d quant_residuals_omega;
};

/// \brief Structure to hold a generic ECEF position measurement and its covariance.
struct PosMeasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the measurement.
    double time;
    /// Position vector in ECEF frame.
    Eigen::Vector3d r_eb_e;
    /// Covariance matrix of the position measurement.
    Eigen::Matrix3d cov_mat;
};

/// \brief Structure to hold a generic ECEF pose (position and rotation) measurement and its covariance.
/// The rotation part in the covariance matrix is for Euler angles.
struct PosRotMeasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the measurement.
    double time;
    /// Position vector in ECEF frame.
    Eigen::Vector3d r_eb_e;
    /// Body-to-ECEF rotation matrix.
    Eigen::Matrix3d C_b_e;
    /// Covariance matrix of the position and rotation measurement.
    Eigen::Matrix<double,6,6> cov_mat;
};

/// \brief Structure to hold Loosely Coupled (LC) GNSS position and velocity measurements.
struct GnssPosVelMeasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the measurement.
    double time;
    /// Position vector of antenna in ECEF frame.
    Eigen::Vector3d r_ea_e; 
    /// Velocity vector of antenna in ECEF frame.
    Eigen::Vector3d v_ea_e;
    /// Covariance matrix of the position and velocity measurement.
    Eigen::Matrix<double,6,6> cov_mat;
};

/// \brief Structure to define IMU error parameters.
struct ImuErrors {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Accelerometer biases (m/s^2).
    Eigen::Vector3d b_a;
    /// Gyro biases (rad/s).
    Eigen::Vector3d b_g;
    /// Accelerometer scale factor and cross coupling errors.
    Eigen::Matrix3d M_a;
    /// Gyro scale factor and cross coupling errors.
    Eigen::Matrix3d M_g;
    /// Gyro g-dependent biases (rad-sec/m).
    Eigen::Matrix3d G_g;
    /// Accelerometer noise root PSD (m s^-1.5).
    double accel_noise_root_psd;
    /// Gyro noise root PSD (rad s^-0.5).
    double gyro_noise_root_psd;
    /// Accelerometer quantization level (m/s^2).
    double accel_quant_level;
    /// Gyro quantization level (rad/s).
    double gyro_quant_level;
};

/// \brief Structure to hold navigation solution in the NED (North, East, Down) frame.
struct NavSolutionNed {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the navigation solution.
    double time;
    /// Latitude in radians.
    double latitude;
    /// Longitude in radians.
    double longitude;
    /// Height in meters.
    double height;
    /// Velocity in NED frame (North, East, Down).
    Eigen::Vector3d v_eb_n;
    /// Body-to-NED rotation matrix.
    Eigen::Matrix3d C_b_n;
};

/// \brief Structure to hold navigation solution in the ECEF (Earth-Centered, Earth-Fixed) frame.
struct NavSolutionEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the navigation solution.
    double time;
    /// Position vector in ECEF frame.
    Eigen::Vector3d r_eb_e;
    /// Velocity vector in ECEF frame.
    Eigen::Vector3d v_eb_e;
    /// Body-to-ECEF rotation matrix.
    Eigen::Matrix3d C_b_e;
};

/// Structure to hold the error state vector for an ECEF navigation filter.
struct StateEstEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Wether update passed checks or not
    bool valid;
    /// The estimated navigation solution in ECEF frame.
    NavSolutionEcef nav_sol;
    /// Estimated accelerometer biases.
    Eigen::Vector3d acc_bias;
    /// Estimated gyroscope biases.
    Eigen::Vector3d gyro_bias;
    // Estimated receiver clock bias.
    double clock_offset;
    /// Estimated receiver clock drift.
    double clock_drift;
    /// The error covariance matrix.
    Eigen::Matrix<double,17,17> P_matrix;
};

/// \brief Structure to hold navigation errors in the NED (North, East, Down) frame.
struct ErrorsNed {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the error calculation.
    double time;
    /// NED Position error.
    Eigen::Vector3d delta_r_eb_n;
    /// NED Velocity error.
    Eigen::Vector3d delta_v_eb_n;
    /// NED orientation error (Euler angles).
    Eigen::Vector3d delta_rot_nb_n;
};

/// \brief Structure to hold navigation errors in ECEF frame along with their estimated standard deviations.
struct ErrorsSigmasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Time of the error calculation.
    double time;
    /// Position error in ECEF frame.
    Eigen::Vector3d delta_r_eb_e;
    /// Velocity error in ECEF frame.
    Eigen::Vector3d delta_v_eb_e;
    /// Orientation error (rotation vector) in ECEF frame.
    Eigen::Vector3d delta_rot_eb_e;
    /// Standard deviation of position error in ECEF frame.
    Eigen::Vector3d sigma_delta_r_eb_e;
    /// Standard deviation of velocity error in ECEF frame.
    Eigen::Vector3d sigma_delta_v_eb_e;
    /// Standard deviation of orientation error (rotation ector) in ECEF frame.
    Eigen::Vector3d sigma_delta_rot_eb_e;
};

/// \brief Structure to configure Kalman Filter parameters.
struct KfConfig {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Initial attitude uncertainty per axis (deg, converted to rad)
    double init_att_unc;
    /// Initial velocity uncertainty per axis (m/s)
    double init_vel_unc;
    /// Initial position uncertainty per axis (m)
    double init_pos_unc;
    /// Initial accelerometer bias uncertainty per instrument (m/s^2)
    double init_b_a_unc;
    /// Initial gyro bias uncertainty per instrument (rad/sec)
    double init_b_g_unc;
    /// Initial clock offset uncertainty per axis (m)
    double init_clock_offset_unc;
    /// Initial clock drift uncertainty per axis (m/s)
    double init_clock_drift_unc;
    /// Gyro noise PSD (rad^2/s)                
    double gyro_noise_psd;
    /// Accelerometer noise PSD (m^2 s^-3)                
    double accel_noise_psd;
    /// Accelerometer bias random walk PSD (m^2 s^-5)
    double accel_bias_psd;
    /// Gyro bias random walk PSD (rad^2 s^-3)
    double gyro_bias_psd;
    /// Receiver clock frequency-drift PSD (m^2/s^3)
    double clock_freq_psd;
    /// Receiver clock phase-drift PSD (m^2/s)
    double clock_phase_psd;
};

/// \brief Structure to configure GNSS simulation parameters.
struct GnssConfig {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Interval between GNSS epochs (s)
    double epoch_interval;
    /// Initial estimated position (m; ECEF)
    Eigen::Vector3d init_est_r_ea_e;
    /// Number of satellites in constellation.
    int no_sat;
    /// Orbital radius of satellites (m)
    double r_os;
    /// Inclination angle of satellites (deg)
    double inclination;
    /// Longitude offset of constellation (deg)
    double const_delta_lambda;
    /// Timing offset of constellation (s)
    double const_delta_t;
    /// Mask angle (deg)
    double mask_angle;
    /// Signal in space error SD (m) *Give residual where corrections are applied
    double sis_err_sd;
    /// Zenith ionosphere error SD (m) *Give residual where corrections are applied
    double zenith_iono_err_sd;
    /// Zenith troposphere error SD (m) *Give residual where corrections are applied
    double zenith_trop_err_sd;
    /// Code tracking error SD (m) *Can extend to account for multipath
    double code_track_err_sd;
    /// Range rate tracking error SD (m/s) *Can extend to account for multipath
    double rate_track_err_sd;
    /// Receiver clock offset at time=0 (m);
    double rx_clock_offset;
    /// Receiver clock drift at time=0 (m/s);
    double rx_clock_drift;
    /// Pos SD for lc integration
    double lc_pos_sd;
    /// Pos SD for lc integration
    double lc_vel_sd;
    /// Pseudo range SD for tight integration
    double pseudo_range_sd;
    /// Pseudo range rate SD for tight integration
    double range_rate_sd;
};

/// \brief Structure to hold GNSS satellite positions and velocities.
struct SatPosVel {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Satellite ECEF positions.
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, kMaxGnssSatellites, 3> sat_r_es_e;
    /// Satellite ECEF velocities.
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, kMaxGnssSatellites, 3> sat_v_es_e;
};

/// \brief Structure to hold GNSS measurement data.
// gnss_measurements     GNSS measurement data:
//     Column 1              Pseudo-range measurements (m)
//     Column 2              Pseudo-range rate measurements (m/s)
//     Columns 3-5           Satellite ECEF position (m)
//     Columns 6-8           Satellite ECEF velocity (m/s)
struct GnssMeasurements {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Matrix containing pseudo-range, pseudo-range rate, satellite position, and satellite velocity.
    Eigen::Matrix<double, Eigen::Dynamic, 8, 0, kMaxGnssSatellites, 8> meas; 
    /// Number of measurements.
    int no_meas;
    /// Covariance matrix of the measurements.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* kMaxGnssSatellites, 2* kMaxGnssSatellites> cov_mat;
};

// Processed GNSS measurements for LC integration (pos + vel + clock offset)
// Result of LS solver
/// \brief Structure to hold processed GNSS measurements for Loosely Coupled (LC) integration.
/// This is typically the result of a Least Squares (LS) solver.
struct GnssLsPosVelClock {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Estimated antenna position in ECEF frame.
    Eigen::Vector3d r_ea_e; 
    /// Estimated antenna velocity in ECEF frame.
    Eigen::Vector3d v_ea_e;
    /// Estimated receiver clock offset and drift.
    Eigen::Vector2d clock;
};

/// @}

};

#endif