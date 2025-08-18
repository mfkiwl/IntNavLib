
#ifndef CONSTANTS_TYPES_H
#define CONSTANTS_TYPES_H

#include <eigen3/Eigen/Dense>

namespace intnavlib {

/// \defgroup constants_types Constants and Types
/// @{

/// Max number of GNSS satellites.
/// Leaving some extra room for future expansion.
constexpr int MAX_GNSS_SATELLITES = 35;

/// Epsilon for checking if a value is approximately zero.
constexpr double EPSILON = 1.0e-100; 

/// WGS84 Equatorial radius in meters.
constexpr double R_0 = 6378137.0;  // WGS84 Equatorial radius in meters
/// WGS84 eccentricity.
constexpr double e = 0.0818191908425; // WGS84 eccentricity
/// Earth rotation rate.
constexpr double omega_ie = 7.292115e-5;  // Earth rotation rate
/// WGS84 Earth gravitational constant (m^3 s^-2).
constexpr double mu = 3.986004418e14; // WGS84 Earth gravitational constant (m^3 s^-2)
/// WGS84 Earth's second gravitational constant.
constexpr double J_2 = 1.082627e-3; // WGS84 Earth's second gravitational constant
/// Speed of light in m/s.
constexpr double c = 299792458.0 ; // Speed of light in m/s

/// Conversion factor from degrees to radians.
constexpr double deg_to_rad = 0.01745329252;
/// Conversion factor from radians to degrees.
constexpr double rad_to_deg = 1.0/deg_to_rad;
/// Conversion factor from micro-g to meters per second squared.
constexpr double micro_g_to_meters_per_second_squared = 9.80665e-6;

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
    double accel_noise_root_PSD;
    /// Gyro noise root PSD (rad s^-0.5).
    double gyro_noise_root_PSD;
    /// Accelerometer quantization level (m/s^2).
    double accel_quant_level;
    /// Gyro quantization level (rad/s).
    double gyro_quant_level;
    
    /// Default constructor: tactical grade IMU
    ImuErrors() {
        b_a << 900.0,-1300.0,800.0;
        b_a = b_a * micro_g_to_meters_per_second_squared;
        b_g << -9.0, 13.0, -8.0;
        b_g = b_g * deg_to_rad / 3600.0;
        M_a << 500.0, -300.0, 200.0,
                        -150.0, -600.0, 250.0,
                        -250.0,  100.0, 450.0;
        M_a = M_a * 1.0e-6;
        M_g << 400.0, -300.0,  250.0,
                            0.0, -300.0, -150.0,
                            0.0,    0.0, -350.0; 
        M_g = M_g * 1.0e-6;
        G_g << 0.9, -1.1, -0.6,
                        -0.5,  1.9, -1.6,
                        0.3,  1.1, -1.3;
        G_g = G_g * deg_to_rad / (3600.0 * 9.80665);  
        accel_noise_root_PSD = 100.0 * micro_g_to_meters_per_second_squared;
        gyro_noise_root_PSD = 0.01 * deg_to_rad / 60.0;
        accel_quant_level = 1.0e-2;
        gyro_quant_level = 2.0e-4;
    }
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

/// \brief Structure to hold the state estimation after a Loosely Coupled (LC) Kalman Filter update.
struct StateEstEcefLc {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Wether update passed checks or not
    bool valid;
    /// The estimated navigation solution in ECEF frame.
    NavSolutionEcef nav_sol;
    /// Estimated accelerometer biases.
    Eigen::Vector3d acc_bias;
    /// Estimated gyroscope biases.
    Eigen::Vector3d gyro_bias;
    /// The error covariance matrix.
    Eigen::Matrix<double,15,15> P_matrix;
};

/// \brief Structure to hold the state estimation after a Tightly Coupled (TC) Kalman Filter update.
struct StateEstEcefTc {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Wether update passed checks or not
    bool valid;
    /// The estimated navigation solution in ECEF frame.
    NavSolutionEcef nav_sol;
    /// Estimated accelerometer biases.
    Eigen::Vector3d acc_bias;
    /// Estimated gyroscope biases.
    Eigen::Vector3d gyro_bias;
    /// Estimated receiver clock offset.
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
    Eigen::Vector3d delta_eul_nb_n;
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
    /// Orientation error (Euler angles) in ECEF frame.
    Eigen::Vector3d delta_eul_eb_e;
    /// Standard deviation of position error in ECEF frame.
    Eigen::Vector3d sigma_delta_r_eb_e;
    /// Standard deviation of velocity error in ECEF frame.
    Eigen::Vector3d sigma_delta_v_eb_e;
    /// Standard deviation of orientation error (Euler angles) in ECEF frame.
    Eigen::Vector3d sigma_delta_eul_eb_e;
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
    double gyro_noise_PSD;
    /// Accelerometer noise PSD (m^2 s^-3)                
    double accel_noise_PSD;
    /// Accelerometer bias random walk PSD (m^2 s^-5)
    double accel_bias_PSD;
    /// Gyro bias random walk PSD (rad^2 s^-3)
    double gyro_bias_PSD;

    /// Receiver clock frequency-drift PSD (m^2/s^3)
    double clock_freq_PSD;
    /// Receiver clock phase-drift PSD (m^2/s)
    double clock_phase_PSD;

    // Default constructor: tuned for tactical grade IMU
    KfConfig() {
        init_att_unc = deg_to_rad * 1.0;
        init_vel_unc = 0.1;
        init_pos_unc = 10.0;
        init_b_a_unc = 1000.0 * micro_g_to_meters_per_second_squared;
        init_b_g_unc = 10.0 * deg_to_rad / 3600.0;
        init_clock_offset_unc = 10.0;
        init_clock_drift_unc = 0.1;
        gyro_noise_PSD = pow(0.02 * deg_to_rad / 60.0, 2.0);
        accel_noise_PSD = pow(200.0 * micro_g_to_meters_per_second_squared, 2.0);
        accel_bias_PSD = 1.0E-7;
        gyro_bias_PSD = 2.0E-12;
        clock_freq_PSD = 1;
        clock_phase_PSD = 1;
    }
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
    double SIS_err_SD;
    /// Zenith ionosphere error SD (m) *Give residual where corrections are applied
    double zenith_iono_err_SD;
    /// Zenith troposphere error SD (m) *Give residual where corrections are applied
    double zenith_trop_err_SD;
    /// Code tracking error SD (m) *Can extend to account for multipath
    double code_track_err_SD;
    /// Range rate tracking error SD (m/s) *Can extend to account for multipath
    double rate_track_err_SD;
    /// Receiver clock offset at time=0 (m);
    double rx_clock_offset;
    /// Receiver clock drift at time=0 (m/s);
    double rx_clock_drift;

    /// SDs for lc integration
    double lc_pos_sd;
    double lc_vel_sd;

    /// SDs for tc integration
    double pseudo_range_sd;
    double range_rate_sd;

    /// Default constructor
    GnssConfig() {
        epoch_interval = 0.5;
        init_est_r_ea_e = Eigen::Vector3d::Zero();
        no_sat = 30.0;
        r_os = 2.656175E7;
        inclination = 55.0;
        const_delta_lambda = 0.0;
        const_delta_t = 0.0;
        mask_angle = 10.0;
        SIS_err_SD = 1.0;
        zenith_iono_err_SD = 2.0;
        zenith_trop_err_SD = 0.2;
        code_track_err_SD = 1.0;
        rate_track_err_SD = 0.02;
        rx_clock_offset = 10000.0;
        rx_clock_drift = 100.0;
        lc_pos_sd = 2.5;
        lc_vel_sd = 0.1;
        pseudo_range_sd = 2.5;
        range_rate_sd = 0.1;
    }
};

/// \brief Structure to hold GNSS satellite positions and velocities.
struct SatPosVel {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    /// Satellite ECEF positions.
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, MAX_GNSS_SATELLITES, 3> sat_r_es_e;
    /// Satellite ECEF velocities.
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, MAX_GNSS_SATELLITES, 3> sat_v_es_e;
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
    Eigen::Matrix<double, Eigen::Dynamic, 8, 0, MAX_GNSS_SATELLITES, 8> meas; 
    /// Number of measurements.
    int no_meas;
    /// Covariance matrix of the measurements.
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* MAX_GNSS_SATELLITES, 2* MAX_GNSS_SATELLITES> cov_mat;
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

/// \brief Structure to configure a generic position sensor.
struct GenericPosSensorConfig {
    /// Standard deviation on position for all axes (m)
    double std_pos;
    /// Interval between pos sensor epochs (s)
    double epoch_interval;
};

/// \brief Structure to configure a generic position and rotation sensor.
struct GenericPosRotSensorConfig {
    /// Standard deviation on position for all axes (m)
    double std_pos;
    /// Standard deviation on rotation for all axes (rad)
    double std_rot;
    /// Interval between pos rot sensor epochs (s)
    double epoch_interval;
};

/// @}

};

#endif