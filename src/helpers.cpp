#include "helpers.h"

namespace helpers {

NavSolutionEcef nedToEcef(const NavSolutionNed& nav_sol_ned) {
    // Extract the inputs
    double L_b = nav_sol_ned.latitude;
    double lambda_b = nav_sol_ned.longitude;
    double h_b = nav_sol_ned.height;
    Eigen::Vector3d v_eb_n = nav_sol_ned.velocity_ned;
    Eigen::Matrix3d C_b_n = nav_sol_ned.C_b_n;

    // Compute transverse radius of curvature
    double sin_lat = std::sin(L_b);
    double cos_lat = std::cos(L_b);
    double sin_long = std::sin(lambda_b);
    double cos_long = std::cos(lambda_b);
    double R_E = R_0 / std::sqrt(1 - (e * sin_lat) * (e * sin_lat));

    // Convert position from curvilinear to Cartesian ECEF coordinates
    Eigen::Vector3d r_eb_e;
    r_eb_e(0) = (R_E + h_b) * cos_lat * cos_long;
    r_eb_e(1) = (R_E + h_b) * cos_lat * sin_long;
    r_eb_e(2) = ((1 - e * e) * R_E + h_b) * sin_lat;

    // Compute the ECEF to NED coordinate transformation matrix
    Eigen::Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
             -sin_long,            cos_long,          0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;

    // Transform velocity from NED to ECEF frame
    Eigen::Vector3d v_eb_e = C_e_n.transpose() * v_eb_n;

    // Transform attitude from NED to ECEF frame
    Eigen::Matrix3d C_b_e = C_e_n.transpose() * C_b_n;

    // Construct the output structure
    NavSolutionEcef nav_sol_ecef;
    nav_sol_ecef.position_ecef = r_eb_e;
    nav_sol_ecef.velocity_ecef = v_eb_e;
    nav_sol_ecef.C_b_e = C_b_e;

    return nav_sol_ecef;
}


Eigen::Matrix3d rpyToR(const double & roll, 
                        const double & pitch, 
                        const double & yaw) {

    // Precompute sines and cosines of Euler angles
    double sin_phi = std::sin(roll);
    double cos_phi = std::cos(roll);
    double sin_theta = std::sin(pitch);
    double cos_theta = std::cos(pitch);
    double sin_psi = std::sin(yaw);
    double cos_psi = std::cos(yaw);

    // Calculate the coordinate transformation matrix using the provided equations
    Eigen::Matrix3d C;
    C(0,0) = cos_theta * cos_psi;
    C(0,1) = cos_theta * sin_psi;
    C(0,2) = -sin_theta;
    C(1,0) = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
    C(1,1) = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
    C(1,2) = sin_phi * cos_theta;
    C(2,0) = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
    C(2,1) = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
    C(2,2) = cos_phi * cos_theta;

    return C;
}

};