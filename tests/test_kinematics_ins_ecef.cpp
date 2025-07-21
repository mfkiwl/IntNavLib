#include <gtest/gtest.h>
#include "intnavlib.h"
#include <random>
#include <string>

using namespace intnavlib;

namespace {

constexpr double kMaxAllowedAvgPositionError = 0.001; // meters
constexpr char test_profile_path[] = "../data/Profile_3.csv";  // Provide your test profile path

// Helper function to run the INS without IMU/init errors and return average position error
double runIdealInsSimulation(const std::string& motion_profile_path) {
    MotionProfileReader reader(motion_profile_path);

    NavSolutionNed true_nav_ned, true_nav_ned_old;
    NavSolutionEcef true_nav_ecef, true_nav_ecef_old;
    NavSolutionEcef est_nav_ecef, est_nav_ecef_old;
    NavSolutionNed est_nav_ned;

    ImuMeasurements true_imu_meas;
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();

    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);
    est_nav_ecef_old = true_nav_ecef_old;

    std::vector<double> position_errors;
    std::mt19937 gen(42); // fixed seed for repeatability

    while (reader.readNextRow(true_nav_ned)) {
        double tor_i = true_nav_ned.time - true_nav_ned_old.time;

        true_nav_ecef = nedToEcef(true_nav_ned);
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // Ideal IMU: just pass true measurements
        ImuMeasurements imu_meas = true_imu_meas;

        est_nav_ecef = navEquationsEcef(est_nav_ecef_old, imu_meas, tor_i);
        est_nav_ned = ecefToNed(est_nav_ecef);

        ErrorsNed errors = calculateErrorsNed(true_nav_ned, est_nav_ned);
        double pos_error = errors.delta_r_eb_n.norm();
        position_errors.push_back(pos_error);

        // Advance states
        true_nav_ecef_old = true_nav_ecef;
        est_nav_ecef_old = est_nav_ecef;
        true_nav_ned_old = true_nav_ned;
        imu_meas_old = imu_meas;
    }

    double sum = 0.0;
    for (double e : position_errors) {
        sum += e;
    }

    return sum / position_errors.size();
}

} // namespace

TEST(InertialNavigation, NoIMUErrorsHasLowPositionError) {

    double avg_position_error = runIdealInsSimulation(test_profile_path);

    EXPECT_LT(avg_position_error, kMaxAllowedAvgPositionError)
        << "Avg position error with ideal IMU should be under "
        << kMaxAllowedAvgPositionError << " m, but was "
        << avg_position_error << " m.";
}
