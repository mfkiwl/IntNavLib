#include <gtest/gtest.h>
#include "intnavlib.h"
#include <random>
#include <string>

using namespace intnavlib;

using ScalarType = double;

using Vector3 = Eigen::Matrix<ScalarType,3,1>;
using Vector2 = Eigen::Matrix<ScalarType,2,1>;
using Matrix3 = Eigen::Matrix<ScalarType,3,3>;

// Typedefs for convenience
using Vector3 = Eigen::Matrix<ScalarType,3,1>;
using Vector2 = Eigen::Matrix<ScalarType,2,1>;
using Matrix3 = Eigen::Matrix<ScalarType,3,3>;
using NavSolutionNed = Types<ScalarType>::NavSolutionNed;
using ImuErrors = Types<ScalarType>::ImuErrors;
using GnssConfig = Types<ScalarType>::GnssConfig;
using KfConfig = Types<ScalarType>::KfConfig;
using FileWriter = Helpers<ScalarType>::FileWriter;
using NavSolutionEcef = Types<ScalarType>::NavSolutionEcef;
using ImuMeasurements = Types<ScalarType>::ImuMeasurements;
using SatPosVel = Types<ScalarType>::SatPosVel;
using GnssMeasurements = Types<ScalarType>::GnssMeasurements;
using NavKF = Navigation<ScalarType>::NavKF;
using PosMeasEcef = Types<ScalarType>::PosMeasEcef;
using StateEstEcef = Types<ScalarType>::StateEstEcef;
using EvalDataEcef = Types<ScalarType>::EvalDataEcef;
using MotionProfileReader = Helpers<ScalarType>::MotionProfileReader;

constexpr char test_profile_path[] = "../apps/nav_sim/motion_profiles/Profile_3.csv";  // Provide your test profile path

// Helper function to run the INS without IMU/init errors and return average position error
ScalarType runIdealInsSimulation(const std::string& motion_profile_path) {
    MotionProfileReader reader(motion_profile_path);

    NavSolutionNed true_nav_ned, true_nav_ned_old;
    NavSolutionEcef true_nav_ecef, true_nav_ecef_old;
    NavSolutionEcef est_nav_ecef, est_nav_ecef_old;

    ImuMeasurements true_imu_meas;

    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);
    est_nav_ecef_old = true_nav_ecef_old;

    double pos_error_sum = 0.0;
    long unsigned int count = 0;
    while (reader.readNextRow(true_nav_ned)) {
        ScalarType tor_i = true_nav_ned.time - true_nav_ned_old.time;

        true_nav_ecef = nedToEcef(true_nav_ned);
        true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // Ideal IMU: just pass true measurements
        ImuMeasurements imu_meas = true_imu_meas;

        est_nav_ecef = navEquationsEcef(est_nav_ecef_old, imu_meas, tor_i);

        ScalarType pos_error = (est_nav_ecef.r_eb_e - true_nav_ecef.r_eb_e).norm();
        pos_error_sum += double(pos_error);
        count++;

        // Advance states
        true_nav_ecef_old = true_nav_ecef;
        est_nav_ecef_old = est_nav_ecef;
        true_nav_ned_old = true_nav_ned;
    }

    return pos_error_sum / double(count);
}


TEST(inertial_navigation, test_ins_ideal) {

    ScalarType avg_position_error = runIdealInsSimulation(test_profile_path);

    // Absolute tolerance depends on type for INS
    double max_pos_error;
    max_pos_error = 0.001; // meters

    EXPECT_LT(avg_position_error, max_pos_error)
        << "Mean absolute position error should be under "
        << max_pos_error << " m, but was "
        << avg_position_error << " m.";

    // always print summary for visibility
    std::cout << "[INFO] Avg position error: "
            << avg_position_error << " m" << std::endl;
}
