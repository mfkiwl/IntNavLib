#include <thread>
#include <chrono>
#include <random>
#include <filesystem>

#include <glog/logging.h>


#include "intnavlib.h"

/// @example ins_ecef.cpp
/// Inertial navigation example in ECEF frame

using namespace intnavlib;

int main(int argc, char** argv)
{   
    google::InitGoogleLogging(argv[0]);
    FLAGS_stderrthreshold = 0; // INFO: 0, WARNING: 1, ERROR: 2, FATAL: 3
    FLAGS_colorlogtostderr = 1;

    if(argc != 2) {
        LOG(ERROR) << "Pass the motion profile path";
        return 1;
    }

    std::string motion_profile_filename_in(argv[1]);
    
    std::string new_directory = "../results";
    if (!std::filesystem::exists(new_directory)) {
        std::filesystem::create_directory(new_directory);
    }
    
    std::string base_filename = std::filesystem::path(motion_profile_filename_in).filename().string();
    std::string filename_without_extension = base_filename.substr(0, base_filename.find_last_of('.'));

    std::string errors_sigmas_filename_out = new_directory + "/" + 
                                    filename_without_extension + "_errors_sigmas.csv";

    // ============== Random gen ==============

    std::random_device rd;
    // Mersenne twister PRNG, initialized with random seed
    std::mt19937 gen(rd());

    // ============== IMU parameters ==============

    // Default constructor: tactical grade IMU
    ImuErrors imu_errors;

    // ============== KF config, used here only for P prop ==============

    KfConfig kf_config;

    // ============== Init profile reader + writers ============

    MotionProfileReader reader(motion_profile_filename_in);
    ErrorsSigmasEcefWriter errors_sigmas_writer(errors_sigmas_filename_out);

    // ============== Simulation state ==============

    // true nav sol
    NavSolutionNed true_nav_ned_t0;
    reader.readNextRow(true_nav_ned_t0);
    NavSolutionEcef true_nav_ecef_old = nedToEcef(true_nav_ned_t0);

    // est nav sol: perturb true nav sol
    NavSolutionEcef est_nav_ecef_old = true_nav_ecef_old;
    std::normal_distribution att_d{0.0, kf_config.init_att_unc};
    std::normal_distribution vel_d{0.0, kf_config.init_vel_unc};
    std::normal_distribution pos_d{0.0, kf_config.init_pos_unc};
    est_nav_ecef_old.C_b_e = true_nav_ecef_old.C_b_e * rpyToR(Eigen::Vector3d(att_d(gen), att_d(gen), att_d(gen)));
    est_nav_ecef_old.r_eb_e += Eigen::Vector3d(pos_d(gen), pos_d(gen), pos_d(gen));
    est_nav_ecef_old.v_eb_e += Eigen::Vector3d(vel_d(gen), vel_d(gen), vel_d(gen));

    //  imu meas
    ImuMeasurements imu_meas_old;
    imu_meas_old.quant_residuals_f = Eigen::Vector3d::Zero();
    imu_meas_old.quant_residuals_omega = Eigen::Vector3d::Zero();

    // Init error state covariance matrix
    Eigen::Matrix<double,15,15> P_matrix = initializeLcPMmatrix(kf_config);

    while (true) {

        // ========= Get ground truth ============

        // Break if done reading input profile
        NavSolutionNed true_nav_ned;
        if(!reader.readNextRow(true_nav_ned)) break;

        double tor_i = true_nav_ned.time - est_nav_ecef_old.time;
        NavSolutionEcef true_nav_ecef = nedToEcef(true_nav_ned);

        // ========== IMU Simulation ==========

        // Get true specific force and angular rates
        ImuMeasurements ideal_imu_meas = kinematicsEcef(true_nav_ecef, true_nav_ecef_old);
        
        // Get imu measurements by applying IMU model
        ImuMeasurements imu_meas = imuModel(ideal_imu_meas, imu_meas_old, imu_errors, tor_i, gen);

        // ========== Nav equations ==========

        NavSolutionEcef est_nav_ecef = navEquationsEcef(est_nav_ecef_old, imu_meas, tor_i);
        NavSolutionNed est_nav_ned = ecefToNed(est_nav_ecef);

        // ========== Prop uncertainty ==========

        P_matrix  = lcPropUnc(P_matrix, 
                            est_nav_ecef,
                            est_nav_ned,
                            imu_meas,
                            kf_config,
                            tor_i);

        // ========== Write Results ==========

        ErrorsSigmasEcef errors_sigmas_ecef; 
        errors_sigmas_ecef.time = true_nav_ned.time;
        errors_sigmas_ecef.delta_r_eb_e = est_nav_ecef.r_eb_e - true_nav_ecef.r_eb_e;
        errors_sigmas_ecef.delta_v_eb_e = est_nav_ecef.v_eb_e - true_nav_ecef.v_eb_e;
        errors_sigmas_ecef.delta_rot_eb_e = deSkew(est_nav_ecef.C_b_e * true_nav_ecef.C_b_e.transpose() - Eigen::Matrix3d::Identity());
        
        errors_sigmas_ecef.sigma_delta_r_eb_e << sqrt(P_matrix(6,6)), sqrt(P_matrix(7,7)), sqrt(P_matrix(8,8));
        errors_sigmas_ecef.sigma_delta_v_eb_e << sqrt(P_matrix(3,3)), sqrt(P_matrix(4,4)), sqrt(P_matrix(5,5));
        errors_sigmas_ecef.sigma_delta_rot_eb_e << sqrt(P_matrix(0,0)), sqrt(P_matrix(1,1)), sqrt(P_matrix(2,2));
        errors_sigmas_writer.writeNextRow(errors_sigmas_ecef);
        
        // ============= Update simulation state ==============

        true_nav_ecef_old = true_nav_ecef;
        est_nav_ecef_old = est_nav_ecef;
        imu_meas_old = imu_meas;

    }

    return 0;
}