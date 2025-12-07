#ifndef SIMULATION_H
#define SIMULATION_H

#include <eigen3/Eigen/Dense>

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    /// \defgroup simulation Simulation
    /// IMU, GNSS, and generic position / pose sensor behavioral models and utilities.
    /// @{

    template<typename T>
    struct Simulation {

        using Vector3 = Eigen::Matrix<T,3,1>;
        using Vector2 = Eigen::Matrix<T,2,1>;
        using Matrix3 = Eigen::Matrix<T,3,3>;

        using ImuMeasurements = typename Types<T>::ImuMeasurements;
        using PosMeasEcef = typename Types<T>::PosMeasEcef;
        using PosRotMeasEcef = typename Types<T>::PosRotMeasEcef;
        using GnssPosVelMeasEcef = typename Types<T>::GnssPosVelMeasEcef;
        using ImuErrors = typename Types<T>::ImuErrors;
        using NavSolutionNed = typename Types<T>::NavSolutionNed;
        using NavSolutionEcef = typename Types<T>::NavSolutionEcef;
        using StateEstEcef = typename Types<T>::StateEstEcef;
        using ErrorsNed = typename Types<T>::ErrorsNed;
        using EvalDataEcef = typename Types<T>::EvalDataEcef;
        using KfConfig = typename Types<T>::KfConfig;
        using GnssConfig = typename Types<T>::GnssConfig;
        using SatPosVel = typename Types<T>::SatPosVel;
        using GnssMeasurements = typename Types<T>::GnssMeasurements;
        using GnssLsPosVelClock = typename Types<T>::GnssLsPosVelClock;

        static constexpr auto kMaxGnssSatellites = Constants<T>::kMaxGnssSatellites;
        static constexpr auto kEpsilon = Constants<T>::kEpsilon;
        static constexpr auto kR0 = Constants<T>::kR0;
        static constexpr auto kEccentricity = Constants<T>::kEccentricity;
        static constexpr auto kOmega_ie = Constants<T>::kOmega_ie;
        static constexpr auto kGravConst = Constants<T>::kGravConst;
        static constexpr auto kJ2 = Constants<T>::kJ2;
        static constexpr auto kC = Constants<T>::kC;
        static constexpr auto kDegToRad = Constants<T>::kDegToRad;
        static constexpr auto kRadToDeg = Constants<T>::kRadToDeg;
        static constexpr auto kMuGToMetersPerSecondSquared = Constants<T>::kMuGToMetersPerSecondSquared;

        /// \brief Calculates true IMU measurements (specific force and angular velocity) from kinematic states.
        /// This function determines the IMU outputs that would be observed given the change in navigation states.
        /// \param[in] new_nav The navigation solution at the current time step.
        /// \param[in] old_nav The navigation solution at the previous time step.
        /// \return An ImuMeasurements structure containing the true specific force and angular velocity.
        static ImuMeasurements kinematicsEcef(const NavSolutionEcef & new_nav,
                                        const NavSolutionEcef & old_nav);

        /// \brief Models IMU sensor outputs by adding errors (biases, noise, quantization) to true IMU measurements.
        /// \param[in] true_imu_meas The true IMU measurements (specific force and angular velocity).
        /// \param[in] old_imu_meas The IMU measurements from the previous time step, used for quantization residuals.
        /// \param[in] imu_errors Structure containing IMU error parameters (biases, scale factors, noise PSDs, quantization levels).
        /// \param[in] tor_i The integration time interval (delta_t) in seconds.
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return An ImuMeasurements structure containing the modeled IMU outputs with errors.
        static ImuMeasurements imuModel(const ImuMeasurements & true_imu_meas, 
                                const ImuMeasurements & old_imu_meas,
                                const ImuErrors & imu_errors,
                                const T & tor_i,
                                std::mt19937 & gen);

        /// \brief Models a generic ECEF position sensor by adding Gaussian noise to the true position.
        /// \param[in] true_nav The true navigation solution in ECEF frame.
        /// \param[in] pos_sigma The standard deviation of the position measurement noise (in meters).
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return A PosMeasEcef structure containing the noisy position measurement and its covariance.
        static PosMeasEcef genericPosSensModel(const NavSolutionEcef & true_nav,
                                    const T & pos_sigma,
                                    std::mt19937 & gen);

        /// \brief Models a generic ECEF pose (position and rotation) sensor by adding Gaussian noise to the true position and rotation.
        /// \param[in] true_nav The true navigation solution in ECEF frame.
        /// \param[in] pos_sigma The standard deviation of the position measurement noise (in meters).
        /// \param[in] rot_sigma The standard deviation of the rotation measurement noise (in radians).
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return A PosRotMeasEcef structure containing the noisy position and rotation measurement and its covariance.
        static PosRotMeasEcef genericPosRotSensModel(const NavSolutionEcef & true_nav,
                                            const T & pos_sigma,
                                            const T & rot_sigma,
                                            std::mt19937 & gen);

        /// \brief Generates satellite positions and velocities in the ECEF frame based on a simplified circular orbit model.
        /// \param[in] time Current simulation time in seconds.
        /// \param[in] gnss_config Configuration parameters for the GNSS constellation.
        /// \return A SatPosVel structure containing the ECEF Cartesian positions and velocities of all satellites.
        static SatPosVel satellitePositionsAndVelocities(const T & time, const GnssConfig & gnss_config);

        /// \brief Generates GNSS pseudo-range and pseudo-range rate measurements for satellites above the elevation mask angle.
        /// This function also includes satellite positions and velocities in the output.
        /// \param[in] time Current simulation time in seconds.
        /// \param[in] gnss_pos_vel Satellite positions and velocities.
        /// \param[in] true_nav_ned True navigation solution in NED frame.
        /// \param[in] true_nav_ecef True navigation solution in ECEF frame.
        /// \param[in] gnss_biases Pre-calculated GNSS range biases for each satellite.
        /// \param[in] gnss_config Configuration parameters for the GNSS receiver and constellation.
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return A GnssMeasurements structure containing the generated pseudo-range and pseudo-range rate measurements,
        ///         along with corresponding satellite positions and velocities.
        static GnssMeasurements generateGnssMeasurements(const T & time,
                                            const SatPosVel & gnss_pos_vel,
                                            const NavSolutionNed& true_nav_ned,
                                            const NavSolutionEcef& true_nav_ecef,
                                            const Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites> & gnss_biases, 
                                            const GnssConfig& gnss_config,
                                            std::mt19937 & gen);

        /// \brief Initializes GNSS range biases due to signal-in-space, ionosphere, and troposphere errors.
        /// These biases are calculated based on satellite elevation angles and configured error standard deviations.
        /// \param[in] true_nav_ecef True navigation solution in ECEF frame.
        /// \param[in] true_nav_ned True navigation solution in NED frame.
        /// \param[in] gnss_pos_vel Satellite positions and velocities.
        /// \param[in] gnss_config Configuration parameters for the GNSS constellation and error models.
        /// \param[in] gen A Mersenne Twister random number generator for noise generation.
        /// \return An Eigen::Matrix containing the initialized GNSS range biases for each satellite.
        Eigen::Matrix<T, Eigen::Dynamic, 1, 0, kMaxGnssSatellites>
        static initializeGnssBiases(const NavSolutionEcef & true_nav_ecef,
                                        const NavSolutionNed & true_nav_ned,
                                        const SatPosVel & gnss_pos_vel,
                                        const GnssConfig& gnss_config,
                                        std::mt19937 & gen);

    };
    
    /// @}

// Convenience wrappers

using Simulationd = Simulation<double>;

using Simulationf = Simulation<float>;

// double

inline Typesd::ImuMeasurements kinematicsEcef(const Typesd::NavSolutionEcef & new_nav,
                                const Typesd::NavSolutionEcef & old_nav) {
    return Simulationd::kinematicsEcef(new_nav, old_nav);
}

inline Typesd::ImuMeasurements imuModel(const Typesd::ImuMeasurements & true_imu_meas, 
                        const Typesd::ImuMeasurements & old_imu_meas,
                        const Typesd::ImuErrors & imu_errors,
                        const double & tor_i,
                        std::mt19937 & gen) {
    return Simulationd::imuModel(true_imu_meas, old_imu_meas, imu_errors, tor_i, gen);
}

inline Typesd::PosMeasEcef genericPosSensModel(const Typesd::NavSolutionEcef & true_nav, 
                            const double & pos_sigma,
                            std::mt19937 & gen) {
    return Simulationd::genericPosSensModel(true_nav, pos_sigma, gen);
}

inline Typesd::PosRotMeasEcef genericPosRotSensModel(const Typesd::NavSolutionEcef & true_nav, 
                                    const double & pos_sigma,
                                    const double & rot_sigma,
                                    std::mt19937 & gen) {
    return Simulationd::genericPosRotSensModel(true_nav, pos_sigma, rot_sigma, gen);
}

inline Typesd::SatPosVel satellitePositionsAndVelocities(const double & time, const Typesd::GnssConfig & gnss_config) {
    return Simulationd::satellitePositionsAndVelocities(time, gnss_config);
}

inline Typesd::GnssMeasurements generateGnssMeasurements(const double & time,
                                    const Typesd::SatPosVel & gnss_pos_vel,
                                    const Typesd::NavSolutionNed& true_nav_ned,
                                    const Typesd::NavSolutionEcef& true_nav_ecef,
                                    const Eigen::Matrix<double, Eigen::Dynamic, 1, 0, Constantsd::kMaxGnssSatellites> & gnss_biases, 
                                    const Typesd::GnssConfig& gnss_config,
                                    std::mt19937 & gen) {
    return Simulationd::generateGnssMeasurements(time, gnss_pos_vel, true_nav_ned, true_nav_ecef, gnss_biases, gnss_config, gen);
}

inline Eigen::Matrix<double, Eigen::Dynamic, 1, 0, Constantsd::kMaxGnssSatellites>
initializeGnssBiases(const Typesd::NavSolutionEcef & true_nav_ecef,
                        const Typesd::NavSolutionNed & true_nav_ned,
                        const Typesd::SatPosVel & gnss_pos_vel,
                        const Typesd::GnssConfig& gnss_config,
                        std::mt19937 & gen) {
    return Simulationd::initializeGnssBiases(true_nav_ecef, true_nav_ned, gnss_pos_vel, gnss_config, gen);
}

// float

inline Typesf::ImuMeasurements kinematicsEcef(const Typesf::NavSolutionEcef & new_nav,
                                const Typesf::NavSolutionEcef & old_nav) {
    return Simulationf::kinematicsEcef(new_nav, old_nav);
}

inline Typesf::ImuMeasurements imuModel(const Typesf::ImuMeasurements & true_imu_meas, 
                        const Typesf::ImuMeasurements & old_imu_meas,
                        const Typesf::ImuErrors & imu_errors,
                        const float & tor_i,
                        std::mt19937 & gen) {
    return Simulationf::imuModel(true_imu_meas, old_imu_meas, imu_errors, tor_i, gen);
}

inline Typesf::PosMeasEcef genericPosSensModel(const Typesf::NavSolutionEcef & true_nav, 
                            const float & pos_sigma,
                            std::mt19937 & gen) {
    return Simulationf::genericPosSensModel(true_nav, pos_sigma, gen);
}

inline Typesf::PosRotMeasEcef genericPosRotSensModel(const Typesf::NavSolutionEcef & true_nav, 
                                    const float & pos_sigma,
                                    const float & rot_sigma,
                                    std::mt19937 & gen) {
    return Simulationf::genericPosRotSensModel(true_nav, pos_sigma, rot_sigma, gen);
}

inline Typesf::SatPosVel satellitePositionsAndVelocities(const float & time, const Typesf::GnssConfig & gnss_config) {
    return Simulationf::satellitePositionsAndVelocities(time, gnss_config);
}

inline Typesf::GnssMeasurements generateGnssMeasurements(const float & time,
                                    const Typesf::SatPosVel & gnss_pos_vel,
                                    const Typesf::NavSolutionNed& true_nav_ned,
                                    const Typesf::NavSolutionEcef& true_nav_ecef,
                                    const Eigen::Matrix<float, Eigen::Dynamic, 1, 0, Constantsf::kMaxGnssSatellites> & gnss_biases, 
                                    const Typesf::GnssConfig& gnss_config,
                                    std::mt19937 & gen) {
    return Simulationf::generateGnssMeasurements(time, gnss_pos_vel, true_nav_ned, true_nav_ecef, gnss_biases, gnss_config, gen);
}

inline Eigen::Matrix<float, Eigen::Dynamic, 1, 0, Constantsf::kMaxGnssSatellites>
initializeGnssBiases(const Typesf::NavSolutionEcef & true_nav_ecef,
                        const Typesf::NavSolutionNed & true_nav_ned,
                        const Typesf::SatPosVel & gnss_pos_vel,
                        const Typesf::GnssConfig& gnss_config,
                        std::mt19937 & gen) {
    return Simulationf::initializeGnssBiases(true_nav_ecef, true_nav_ned, gnss_pos_vel, gnss_config, gen);
}

extern template struct Simulation<double>;
extern template struct Simulation<float>;

};

#endif // SIMULATION_H