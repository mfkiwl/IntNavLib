#ifndef HELPERS_H
#define HELPERS_H

#include <random>
#include <ctime>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include "constants_types.h"

namespace intnavlib {
    
    // Gets gravity vector at given ecef position, in ecef frame
    Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e);

    // NED to ECEF
    NavSolutionEcef nedToEcef(const NavSolutionNed & nav_sol_ned);

    // ECEF to NED
    NavSolutionNed ecefToNed(const NavSolutionEcef & nav_sol_ecef);

    // Converts degrees to radians
    inline double degToRad(const double & degrees) {
        return degrees * 0.01745329252;
    }

    // Converts radians to degrees
    inline double radToDeg(const double & rads) {
        return rads / 0.01745329252;
    }

    Eigen::Matrix3d rpyToR(const Eigen::Vector3d & rpy);

    Eigen::Vector3d rToRpy(const Eigen::Matrix3d & C);

    Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d & a);

    // signum function
    template <typename T> inline T sgn(const T & val) {
    return T((T(0) < val) - (val < T(0)));
    }

    // Function to get the current date and time as a formatted string
    std::string getCurrentDateTime();

    // Get radii of curvature from latitude
    Eigen::Vector2d radiiOfCurvature(double L);

    ErrorsNed calculateErrorsNed(const NavSolutionNed & true_nav_sol, 
                                const NavSolutionNed & est_nav_sol);

};

#endif