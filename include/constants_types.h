
#ifndef CONSTANTS_TYPES_H
#define CONSTANTS_TYPES_H

#include <eigen3/Eigen/Dense>

namespace intnavlib {

// Max number of gnss satellites
// Leaving some extra room
constexpr int MAX_GNSS_SATELLITES = 35;

// Epsilon for == 0 checks
// Todo tune
constexpr double EPSILON = 1.0e-100; 

// Constants for WGS84 model
constexpr double R_0 = 6378137.0;  // WGS84 Equatorial radius in meters
constexpr double e = 0.0818191908425; // WGS84 eccentricity
constexpr double omega_ie = 7.292115e-5;  // Earth rotation rate
constexpr double mu = 3.986004418e14; // WGS84 Earth gravitational constant (m^3 s^-2)
constexpr double J_2 = 1.082627e-3; // WGS84 Earth's second gravitational constant
constexpr double c = 299792458.0 ; // Speed of light in m/s

constexpr double deg_to_rad = 0.01745329252;
constexpr double rad_to_deg = 1.0/deg_to_rad;
constexpr double micro_g_to_meters_per_second_squared = 9.80665e-6;

// IMU measurements
struct ImuMeasurements {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    Eigen::Vector3d f; // Specific forces
    Eigen::Vector3d omega; // angular velocities
    Eigen::Vector3d quant_residuals_f;
    Eigen::Vector3d quant_residuals_omega;
};

// Generic ECEF position measurement + covariance
struct PosMeasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    Eigen::Vector3d r_eb_e;
    Eigen::Matrix3d cov_mat;
};

// Generic ECEF pose measurement + covariance
// rotation part in cov is for euler angles
struct PosRotMeasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    Eigen::Vector3d r_eb_e;
    Eigen::Matrix3d C_b_e;
    Eigen::Matrix<double,6,6> cov_mat;
};

// LC GNSS pos + vel meas
struct GnssPosVelMeasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Vector3d r_ea_e; 
    Eigen::Vector3d v_ea_e;
    Eigen::Matrix<double,6,6> cov_mat;
};

// IMU errors
struct ImuErrors {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Vector3d b_a;              // Accelerometer biases (m/s^2)
    Eigen::Vector3d b_g;              // Gyro biases (rad/s)
    Eigen::Matrix3d M_a;              // Accelerometer scale factor and cross coupling errors
    Eigen::Matrix3d M_g;              // Gyro scale factor and cross coupling errors            
    Eigen::Matrix3d G_g;              // Gyro g-dependent biases (rad-sec/m)             
    double accel_noise_root_PSD;   // Accelerometer noise root PSD (m s^-1.5)
    double gyro_noise_root_PSD;    // Gyro noise root PSD (rad s^-0.5)
    double accel_quant_level;      // Accelerometer quantization level (m/s^2)
    double gyro_quant_level;       // Gyro quantization level (rad/s)
};

// Nav solution in NED frame
struct NavSolutionNed {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    double latitude;  // Latitude in radians
    double longitude; // Longitude in radians
    double height;    // Height in meters
    Eigen::Vector3d v_eb_n; // Velocity in NED frame (North, East, Down)
    Eigen::Matrix3d C_b_n; // Body-to-NED rotation matrix
};

// Nav solution in ECEF frame
struct NavSolutionEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    Eigen::Vector3d r_eb_e; // Position in ECEF frame
    Eigen::Vector3d v_eb_e; // Velocity in ECEF frame
    Eigen::Matrix3d C_b_e; // Body-to-ECEF rotation matrix
};

// State estimation after KF update
struct StateEstEcefLc {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    NavSolutionEcef nav_sol;
    Eigen::Vector3d acc_bias;
    Eigen::Vector3d gyro_bias;
    Eigen::Matrix<double,15,15> P_matrix;
};

// State estimation after KF update
struct StateEstEcefTc {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    NavSolutionEcef nav_sol;
    Eigen::Vector3d acc_bias;
    Eigen::Vector3d gyro_bias;
    double clock_offset;
    double clock_drift;
    Eigen::Matrix<double,17,17> P_matrix;
};

struct ErrorsNed {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    Eigen::Vector3d delta_r_eb_n; // NED Position error
    Eigen::Vector3d delta_v_eb_n; // NED Velocity error
    Eigen::Vector3d delta_eul_nb_n; // NED orientation error
};

// Errors in ecef + estimated std devs
struct ErrorsSigmasEcef {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double time;
    Eigen::Vector3d delta_r_eb_e; 
    Eigen::Vector3d delta_v_eb_e; 
    Eigen::Vector3d delta_eul_eb_e;
    Eigen::Vector3d sigma_delta_r_eb_e; 
    Eigen::Vector3d sigma_delta_v_eb_e; 
    Eigen::Vector3d sigma_delta_eul_eb_e;
};

struct KfConfig {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // Initial attitude uncertainty per axis (deg, converted to rad)
    double init_att_unc;
    // Initial velocity uncertainty per axis (m/s)
    double init_vel_unc;
    // Initial position uncertainty per axis (m)
    double init_pos_unc;
    // Initial accelerometer bias uncertainty per instrument (micro-g, converted
    // to m/s^2)
    double init_b_a_unc;
    // Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
    double init_b_g_unc;
    // Initial clock offset uncertainty per axis (m)
    double init_clock_offset_unc;
    // Initial clock drift uncertainty per axis (m/s)
    double init_clock_drift_unc;

    // Gyro noise PSD (deg^2 per hour, converted to rad^2/s)                
    double gyro_noise_PSD;
    // Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)                
    double accel_noise_PSD;
    // Accelerometer bias random walk PSD (m^2 s^-5)
    double accel_bias_PSD;
    // Gyro bias random walk PSD (rad^2 s^-3)
    double gyro_bias_PSD;

    // Receiver clock frequency-drift PSD (m^2/s^3)
    double clock_freq_PSD;
    // Receiver clock phase-drift PSD (m^2/s)
    double clock_phase_PSD;

};

struct GnssConfig {
    // Interval between GNSS epochs (s)
    double epoch_interval;

    // Initial estimated position (m; ECEF)
    Eigen::Vector3d init_est_r_ea_e;

    // Number of satellites in constellation
    int no_sat;
    // Orbital radius of satellites (m)
    double r_os;
    // Inclination angle of satellites (deg)
    double inclination;
    // Longitude offset of constellation (deg)
    double const_delta_lambda;
    // Timing offset of constellation (s)
    double const_delta_t;

    // Mask angle (deg)
    double mask_angle;
    // Signal in space error SD (m) *Give residual where corrections are applied
    double SIS_err_SD;
    // Zenith ionosphere error SD (m) *Give residual where corrections are applied
    double zenith_iono_err_SD;
    // Zenith troposphere error SD (m) *Give residual where corrections are applied
    double zenith_trop_err_SD;
    // Code tracking error SD (m) *Can extend to account for multipath
    double code_track_err_SD;
    // Range rate tracking error SD (m/s) *Can extend to account for multipath
    double rate_track_err_SD;
    // Receiver clock offset at time=0 (m);
    double rx_clock_offset;
    // Receiver clock drift at time=0 (m/s);
    double rx_clock_drift;

    // SDs for lc integration
    double lc_pos_sd;
    double lc_vel_sd;

    // SDs for tc integration
    double pseudo_range_sd;
    double range_rate_sd;
};

// GNSS satellites positions and velocities
struct SatPosVel {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, MAX_GNSS_SATELLITES, 3> sat_r_es_e;
    Eigen::Matrix<double, Eigen::Dynamic, 3, 0, MAX_GNSS_SATELLITES, 3> sat_v_es_e;
};

// gnss_measurements     GNSS measurement data:
//     Column 1              Pseudo-range measurements (m)
//     Column 2              Pseudo-range rate measurements (m/s)
//     Columns 3-5           Satellite ECEF position (m)
//     Columns 6-8           Satellite ECEF velocity (m/s)
struct GnssMeasurements {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Matrix<double, Eigen::Dynamic, 8, 0, MAX_GNSS_SATELLITES, 8> meas; 
    int no_meas;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, 2* MAX_GNSS_SATELLITES, 2* MAX_GNSS_SATELLITES> cov_mat;
};

// Processed GNSS measurements for LC integration (pos + vel + clock offset)
// Result of LS solver
struct GnssLsPosVelClock {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Eigen::Vector3d r_ea_e; 
    Eigen::Vector3d v_ea_e;
    Eigen::Vector2d clock;
};

struct GenericPosSensorConfig {
    // Standard deviation on position for all axes (m)
    double std_pos;
    // Interval between pos sensor epochs (s)
    double epoch_interval;
};

struct GenericPosRotSensorConfig {
    // Standard deviation on position for all axes (m)
    double std_pos;
    // Standard deviation on rotation for all axes (rad)
    double std_rot;
    // Interval between pos rot sensor epochs (s)
    double epoch_interval;
};

}

#endif