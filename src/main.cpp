#include <thread>
#include <chrono>

#include "helpers.h"
#include "lc_ins_gnss.h"
#include "motion_profile_reader.h"

int main(int argc, char** argv)
{   
    MotionProfileReader reader("/home/tommaso/Workspace/grovesCpp/data/Profile_4.csv");

    // Ground truth nav solution in ned
    NavSolutionNed true_nav_ned;
    NavSolutionNed true_nav_ned_old;

    // Ground truth nav solution in ecef
    NavSolutionEcef true_nav_ecef;
    NavSolutionEcef true_nav_ecef_old;

    // Estimated nav solution in ecef
    NavSolutionEcef est_nav_ecef;
    NavSolutionEcef est_nav_ecef_old;

    // Ground truth imu measurements from kinematics
    ImuMeasurements true_imu_meas;

    // Simulated imu measurements
    ImuMeasurements imu_meas;

    // Imu dt
    double dt;

    // Init both true old nav sol, and estimated old nav sol
    reader.readNextRow(true_nav_ned_old);
    true_nav_ecef_old = nedToEcef(true_nav_ned_old);
    est_nav_ecef_old = true_nav_ecef_old; // Should instead use gnss for this

    while (reader.readNextRow(true_nav_ned)) {
        
        // Get true ecef nav sol from motion profile in ned
        true_nav_ecef = nedToEcef(true_nav_ned);

        // // Get true specific force and angular rates
        // true_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);

        // // Get imu measurements by applying IMU model
        // imu_meas = imuModel(true_imu_meas, imu_errors);

        // // Predict ecef nav solution (INS)
        // est_nav_ecef = navEquationsEcef(dt, est_nav_ecef_old, imu_meas);

        // true_nav_ecef_old = true_nav_ecef;
        // est_nav_ecef_old = est_nav_ecef;

    }

    return 0;
}