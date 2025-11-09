
#ifndef CONSTANTS_TYPES_H
#define CONSTANTS_TYPES_H

#include <eigen3/Eigen/Dense>

namespace intnavlib {

/// \defgroup constants_types Constants and Types
/// @{

template <typename T>
struct Constants {
    /// Max number of GNSS satellites.
    /// Leaving some extra room for future expansion.
    static constexpr int kMaxGnssSatellites = 35;
    /// Epsilon for checking if a value is approximately zero.
    static constexpr T kEpsilon = 1.0e-100; 
    /// WGS84 Equatorial radius in meters.
    static constexpr T kR0 = 6378137.0;  // WGS84 Equatorial radius in meters
    /// WGS84 eccentricity.
    static constexpr T kEccentricity = 0.0818191908425; // WGS84 eccentricity
    /// Earth rotation rate.
    static constexpr T kOmega_ie = 7.292115e-5;  // Earth rotation rate
    /// WGS84 Earth gravitational constant (m^3 s^-2).
    static constexpr T kGravConst = 3.986004418e14; // WGS84 Earth gravitational constant (m^3 s^-2)
    /// WGS84 Earth's second gravitational constant.
    static constexpr T kJ2 = 1.082627e-3; // WGS84 Earth's second gravitational constant
    /// Speed of light in m/s.
    static constexpr T kC = 299792458.0 ; // Speed of light in m/s
    /// Conversion factor from degrees to radians.
    static constexpr T kDegToRad = 0.01745329252;
    /// Conversion factor from radians to degrees.
    static constexpr T kRadToDeg = 1.0/kDegToRad;
    /// Conversion factor from micro-g to meters per second squared.
    static constexpr T kMuGToMetersPerSecondSquared = 9.80665e-6;
};

/// \brief Structure to hold IMU measurements.
template <typename T>
struct Types {

    using Vector2 = Eigen::Matrix<T,2,1>;
    using Vector3 = Eigen::Matrix<T,3,1>;
    using Matrix3 = Eigen::Matrix<T,3,3>;

    struct ImuMeasurements {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the measurement.
        T time;
        /// Specific forces (accelerometer measurements).
        Vector3 f;
        /// Angular velocities (gyroscope measurements).
        Vector3 omega;
        /// Quantization residuals for specific forces.
        Vector3 quant_residuals_f;
        /// Quantization residuals for angular velocities.
        Vector3 quant_residuals_omega;
    };

    /// \brief Structure to hold a generic ECEF position measurement and its covariance.
    struct PosMeasEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the measurement.
        T time;
        /// Position vector in ECEF frame.
        Vector3 r_eb_e;
        /// Covariance matrix of the position measurement.
        Matrix3 cov_mat;
    };

    /// \brief Structure to hold a generic ECEF pose (position and rotation) measurement and its covariance.
    /// The rotation part in the covariance matrix is for Euler angles.
    struct PosRotMeasEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the measurement.
        T time;
        /// Position vector in ECEF frame.
        Vector3 r_eb_e;
        /// Body-to-ECEF rotation matrix.
        Matrix3 C_b_e;
        /// Covariance matrix of the position and rotation measurement.
        Eigen::Matrix<T,6,6> cov_mat;
    };

    /// \brief Structure to hold Loosely Coupled (LC) GNSS position and velocity measurements.
    struct GnssPosVelMeasEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the measurement.
        T time;
        /// Position vector of antenna in ECEF frame.
        Vector3 r_ea_e; 
        /// Velocity vector of antenna in ECEF frame.
        Vector3 v_ea_e;
        /// Covariance matrix of the position and velocity measurement.
        Eigen::Matrix<T,6,6> cov_mat;
    };

    /// \brief Structure to define IMU error parameters.
    struct ImuErrors {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Accelerometer biases (m/s^2).
        Vector3 b_a;
        /// Gyro biases (rad/s).
        Vector3 b_g;
        /// Accelerometer scale factor and cross coupling errors.
        Matrix3 M_a;
        /// Gyro scale factor and cross coupling errors.
        Matrix3 M_g;
        /// Gyro g-dependent biases (rad-sec/m).
        Matrix3 G_g;
        /// Accelerometer noise root PSD (m s^-1.5).
        T accel_noise_root_psd;
        /// Gyro noise root PSD (rad s^-0.5).
        T gyro_noise_root_psd;
        /// Accelerometer quantization level (m/s^2).
        T accel_quant_level;
        /// Gyro quantization level (rad/s).
        T gyro_quant_level;
    };

    /// \brief Structure to hold navigation solution in the NED (North, East, Down) frame.
    struct NavSolutionNed {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the navigation solution.
        T time;
        /// Latitude in radians.
        T latitude;
        /// Longitude in radians.
        T longitude;
        /// Height in meters.
        T height;
        /// Velocity in NED frame (North, East, Down).
        Vector3 v_eb_n;
        /// Body-to-NED rotation matrix.
        Matrix3 C_b_n;
    };

    /// \brief Structure to hold navigation solution in the ECEF (Earth-Centered, Earth-Fixed) frame.
    struct NavSolutionEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the navigation solution.
        T time;
        /// Position vector in ECEF frame.
        Vector3 r_eb_e;
        /// Velocity vector in ECEF frame.
        Vector3 v_eb_e;
        /// Body-to-ECEF rotation matrix.
        Matrix3 C_b_e;
    };

    /// Structure to hold the error state vector for an ECEF navigation filter.
    struct StateEstEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Wether update passed checks or not
        bool valid;
        /// The estimated navigation solution in ECEF frame.
        NavSolutionEcef nav_sol;
        /// Estimated accelerometer biases.
        Vector3 acc_bias;
        /// Estimated gyroscope biases.
        Vector3 gyro_bias;
        // Estimated receiver clock bias.
        T clock_offset;
        /// Estimated receiver clock drift.
        T clock_drift;
        /// The error covariance matrix.
        Eigen::Matrix<T,17,17> P_matrix;
        /// Innovations + sigmas vector (dynamic, depends on measurements)
        std::vector<std::pair<T, T>> innovations_sigmas;
        // /// Default constructor
        StateEstEcef()
            : valid(false),
            nav_sol(), // assumes NavSolutionEcef has its own default constructor
            acc_bias(Vector3::Constant(-1.0)),
            gyro_bias(Vector3::Constant(-1.0)),
            clock_offset(-1.0),
            clock_drift(-1.0),
            P_matrix(Eigen::Matrix<T,17,17>::Constant(-1.0)),
            innovations_sigmas()
        {}
    };

    /// \brief Structure to hold navigation errors in the NED (North, East, Down) frame.
    struct ErrorsNed {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the error calculation.
        T time;
        /// NED Position error.
        Vector3 delta_r_eb_n;
        /// NED Velocity error.
        Vector3 delta_v_eb_n;
        /// NED orientation error (Euler angles).
        Vector3 delta_rot_nb_n;
    };

    /// \brief Structure to hold navigation eval data in ECEF frame along with their estimated standard deviations.
    struct EvalDataEcef {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Time of the error calculation.
        T time;
        /// Position error in ECEF frame.
        Vector3 delta_r_eb_e;
        /// Velocity error in ECEF frame.
        Vector3 delta_v_eb_e;
        /// Orientation error (rotation vector) in ECEF frame.
        Vector3 delta_rot_eb_e;
        /// Accelerometer bias estimate
        Vector3 delta_b_a;
        /// Gyro bias estimate
        Vector3 delta_b_g;
        /// Clock bias estimate
        T delta_clock_offset;
        /// Clock drift estimate
        T delta_clock_drift;
        /// Standard deviation of position error in ECEF frame.
        Vector3 sigma_delta_r_eb_e;
        /// Standard deviation of velocity error in ECEF frame.
        Vector3 sigma_delta_v_eb_e;
        /// Standard deviation of orientation error (rotation ector) in ECEF frame.
        Vector3 sigma_delta_rot_eb_e;
        /// Standard deviation of accelerometer bias
        Vector3 sigma_delta_b_a;
        /// Standard deviation of gyro bias
        Vector3 sigma_delta_b_g;
        /// Standard deviation of clock bias
        T sigma_delta_clock_offset;
        /// Standard deviation of clock drift   
        T sigma_delta_clock_drift;
        // Innovations and sigmas
        std::vector<std::pair<T, T>> innovations_sigmas;
    };

    /// \brief Structure to configure Kalman Filter parameters.
    struct KfConfig {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Initial attitude uncertainty per axis (deg, converted to rad)
        T init_att_unc;
        /// Initial velocity uncertainty per axis (m/s)
        T init_vel_unc;
        /// Initial position uncertainty per axis (m)
        T init_pos_unc;
        /// Initial accelerometer bias uncertainty per instrument (m/s^2)
        T init_b_a_unc;
        /// Initial gyro bias uncertainty per instrument (rad/sec)
        T init_b_g_unc;
        /// Initial clock offset uncertainty per axis (m)
        T init_clock_offset_unc;
        /// Initial clock drift uncertainty per axis (m/s)
        T init_clock_drift_unc;
        /// Gyro noise PSD (rad^2/s)                
        T gyro_noise_psd;
        /// Accelerometer noise PSD (m^2 s^-3)                
        T accel_noise_psd;
        /// Accelerometer bias random walk PSD (m^2 s^-5)
        T accel_bias_psd;
        /// Gyro bias random walk PSD (rad^2 s^-3)
        T gyro_bias_psd;
        /// Receiver clock frequency-drift PSD (m^2/s^3)
        T clock_freq_psd;
        /// Receiver clock phase-drift PSD (m^2/s)
        T clock_phase_psd;
    };

    /// \brief Structure to configure GNSS simulation parameters.
    struct GnssConfig {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Interval between GNSS epochs (s)
        T epoch_interval;
        /// Initial estimated position (m; ECEF)
        Vector3 init_est_r_ea_e;
        /// Number of satellites in constellation.
        int no_sat;
        /// Orbital radius of satellites (m)
        T r_os;
        /// Inclination angle of satellites (deg)
        T inclination;
        /// Longitude offset of constellation (deg)
        T const_delta_lambda;
        /// Timing offset of constellation (s)
        T const_delta_t;
        /// Mask angle (deg)
        T mask_angle;
        /// Signal in space error SD (m) *Give residual where corrections are applied
        T sis_err_sd;
        /// Zenith ionosphere error SD (m) *Give residual where corrections are applied
        T zenith_iono_err_sd;
        /// Zenith troposphere error SD (m) *Give residual where corrections are applied
        T zenith_trop_err_sd;
        /// Code tracking error SD (m) *Can extend to account for multipath
        T code_track_err_sd;
        /// Range rate tracking error SD (m/s) *Can extend to account for multipath
        T rate_track_err_sd;
        /// Receiver clock offset at time=0 (m);
        T rx_clock_offset;
        /// Receiver clock drift at time=0 (m/s);
        T rx_clock_drift;
        /// Pos SD for lc integration
        T lc_pos_sd;
        /// Pos SD for lc integration
        T lc_vel_sd;
        /// Pseudo range SD for tight integration
        T pseudo_range_sd;
        /// Pseudo range rate SD for tight integration
        T range_rate_sd;
    };

    /// \brief Structure to hold GNSS satellite positions and velocities.
    struct SatPosVel {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Satellite ECEF positions.
        Eigen::Matrix<T, Eigen::Dynamic, 3, 0, Constants<T>::kMaxGnssSatellites, 3> sat_r_es_e;
        /// Satellite ECEF velocities.
        Eigen::Matrix<T, Eigen::Dynamic, 3, 0, Constants<T>::kMaxGnssSatellites, 3> sat_v_es_e;
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
        Eigen::Matrix<T, Eigen::Dynamic, 8, 0, Constants<T>::kMaxGnssSatellites, 8> meas; 
        /// Number of measurements.
        int no_meas;
        /// Covariance matrix of the measurements.
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, 0, 2* Constants<T>::kMaxGnssSatellites, 2* Constants<T>::kMaxGnssSatellites> cov_mat;
    };

    // Processed GNSS measurements for LC integration (pos + vel + clock offset)
    // Result of LS solver
    /// \brief Structure to hold processed GNSS measurements for Loosely Coupled (LC) integration.
    /// This is typically the result of a Least Squares (LS) solver.
    struct GnssLsPosVelClock {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        /// Estimated antenna position in ECEF frame.
        Vector3 r_ea_e; 
        /// Estimated antenna velocity in ECEF frame.
        Vector3 v_ea_e;
        /// Estimated receiver clock offset and drift.
        Vector2 clock;
    };

};

// For backward compatibility

// Define default type
#ifdef USE_FLOAT
    using nav_type = float;
#else
    using nav_type = double;
#endif

// Constants

using Constantst = Constants<nav_type>;

static constexpr int kMaxGnssSatellites = Constantst::kMaxGnssSatellites;
static constexpr nav_type kEpsilon = Constantst::kEpsilon;
static constexpr nav_type kR0 = Constantst::kR0;
static constexpr nav_type kEccentricity = Constantst::kEccentricity;
static constexpr nav_type kOmega_ie = Constantst::kOmega_ie;
static constexpr nav_type kGravConst = Constantst::kGravConst;
static constexpr nav_type kJ2 = Constantst::kJ2;
static constexpr nav_type kC = Constantst::kC;
static constexpr nav_type kDegToRad = Constantst::kDegToRad;
static constexpr nav_type kRadToDeg = Constantst::kRadToDeg;
static constexpr nav_type kMuGToMetersPerSecondSquared = Constantst::kMuGToMetersPerSecondSquared;

// Types

using Typest = Types<nav_type>;

using ImuMeasurements = Typest::ImuMeasurements;
using PosMeasEcef = Typest::PosMeasEcef;
using PosRotMeasEcef = Typest::PosRotMeasEcef;
using GnssPosVelMeasEcef = Typest::GnssPosVelMeasEcef;
using ImuErrors = Typest::ImuErrors;
using NavSolutionNed = Typest::NavSolutionNed;
using NavSolutionEcef = Typest::NavSolutionEcef;
using StateEstEcef = Typest::StateEstEcef;
using ErrorsNed = Typest::ErrorsNed;
using EvalDataEcef = Typest::EvalDataEcef;
using KfConfig = Typest::KfConfig;
using GnssConfig = Typest::GnssConfig;
using SatPosVel = Typest::SatPosVel;
using GnssMeasurements = Typest::GnssMeasurements;
using GnssLsPosVelClock = Typest::GnssLsPosVelClock;


/// @}

};

#endif // CONSTANTS_TYPES_H