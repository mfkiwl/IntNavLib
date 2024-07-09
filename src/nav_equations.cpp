
#include "nav_equations.h"

namespace intnavlib {

NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i) {

    NavSolutionEcef new_nav;

    new_nav.time = imu_meas.time;

    // ATTITUDE UPDATE
    // From (2.145) determine the Earth rotation over the update interval
    // C_Earth = C_e_i' * old_C_e_i

    double alpha_ie = omega_ie * tor_i;

    double cos_alpha_ie = cos(alpha_ie);
    double sin_alpha_ie = sin(alpha_ie);

    Eigen::Matrix3d C_Earth;
    C_Earth << cos_alpha_ie, sin_alpha_ie, 0.0,
                -sin_alpha_ie, cos_alpha_ie, 0.0,
                            0.0,             0.0,  1.0;
                        
    // Calculate attitude increment, magnitude, and skew-symmetric matrix
    Eigen::Vector3d alpha_ib_b = imu_meas.omega * tor_i;
    double mag_alpha = alpha_ib_b.norm();
    Eigen::Matrix3d Alpha_ib_b = skewSymmetric(alpha_ib_b);  

    // Obtain coordinate transformation matrix from the new attitude w.r.t. an
    // inertial frame to the old using Rodrigues' formula, (5.73)
    Eigen::Matrix3d C_new_old;
    if (mag_alpha>1.0e-8) {
        C_new_old = Eigen::Matrix3d::Identity() + 
                    ((sin(mag_alpha) / mag_alpha) * Alpha_ib_b) +
                    ((1.0 - cos(mag_alpha)) / pow(mag_alpha,2.0)) * Alpha_ib_b * Alpha_ib_b;
    }
    else {
        C_new_old = Eigen::Matrix3d::Identity() + Alpha_ib_b;
    }

    // Update attitude using (5.75)
    new_nav.C_b_e = C_Earth * old_nav.C_b_e * C_new_old;
        
    // SPECIFIC FORCE FRAME TRANSFORMATION
    // Calculate the average body-to-ECEF-frame coordinate transformation
    // matrix over the update interval using (5.84) and (5.85)

    Eigen::Vector3d alpha_ie_vec;
    alpha_ie_vec << 0.0 , 0.0 , alpha_ie;   

    Eigen::Matrix3d ave_C_b_e;
    if (mag_alpha>1.0e-8) {
        ave_C_b_e = old_nav.C_b_e * 
            (Eigen::Matrix3d::Identity() + 

            ((1.0 - cos(mag_alpha)) / pow(mag_alpha,2.0)) *
            Alpha_ib_b + 
            
            ((1.0 - sin(mag_alpha) / mag_alpha) / pow(mag_alpha,2.0)) * 
            Alpha_ib_b * Alpha_ib_b) - 

            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }
    else { // Approximate if angle small enough (sum not multiply)
        ave_C_b_e = old_nav.C_b_e -
            0.5 * skewSymmetric(alpha_ie_vec) * old_nav.C_b_e;
    }

    // Transform specific force to ECEF-frame resolving axes using (5.85)
    Eigen::Vector3d f_ib_e = ave_C_b_e * imu_meas.f;
        
    // UPDATE VELOCITY
    // From (5.36)

    Eigen::Vector3d omega_ie_vec;
    omega_ie_vec << 0.0 , 0.0 , omega_ie;   

    new_nav.v_eb_e = old_nav.v_eb_e + tor_i * (f_ib_e + gravityEcef(old_nav.r_eb_e) -
                            2.0 * skewSymmetric(omega_ie_vec) * old_nav.v_eb_e);


    // UPDATE CARTESIAN POSITION
    // From (5.38),
    new_nav.r_eb_e = old_nav.r_eb_e + (new_nav.v_eb_e + old_nav.v_eb_e) * 0.5 * tor_i; 

    return new_nav;

}

}