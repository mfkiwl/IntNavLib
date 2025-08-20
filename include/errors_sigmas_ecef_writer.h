#ifndef ERRORS_SIGMAS_ECEF_WRITER_H
#define ERRORS_SIGMAS_ECEF_WRITER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "helpers.h"

namespace intnavlib {

class ErrorsSigmasEcefWriter {
public:
    ErrorsSigmasEcefWriter(const std::string& filename) : file(filename), ok(false) {
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        ok = true;
    }

    bool writeNextRow(const ErrorsSigmasEcef & row) {
        if (!file.good()) return false;

         // Extract position errors (in meters)
        Eigen::Vector3d delta_r_eb_e = row.delta_r_eb_e;
        Eigen::Vector3d sigma_delta_r_eb_e = row.sigma_delta_r_eb_e;

        // Extract velocity errors (in m/s)
        Eigen::Vector3d delta_v_eb_e = row.delta_v_eb_e;
        Eigen::Vector3d sigma_delta_v_eb_e = row.sigma_delta_v_eb_e;

        // Convert Euler angles (in radians) to degrees
        Eigen::Vector3d delta_rot_eb_e = row.delta_rot_eb_e * rad_to_deg;
        Eigen::Vector3d sigma_delta_rot_eb_e = row.sigma_delta_rot_eb_e * rad_to_deg;

        // Write the row to the CSV file
        file << row.time << ","

             << delta_r_eb_e[0] << "," // Errors
             << delta_r_eb_e[1] << ","
             << delta_r_eb_e[2] << ","
             << delta_v_eb_e[0] << "," 
             << delta_v_eb_e[1] << ","
             << delta_v_eb_e[2] << ","
             << delta_rot_eb_e[0] << ","
             << delta_rot_eb_e[1] << ","
             << delta_rot_eb_e[2] << ","

             << sigma_delta_r_eb_e[0] << "," // Errors sigma
             << sigma_delta_r_eb_e[1] << ","
             << sigma_delta_r_eb_e[2] << ","
             << sigma_delta_v_eb_e[0] << "," 
             << sigma_delta_v_eb_e[1] << ","
             << sigma_delta_v_eb_e[2] << ","
             << sigma_delta_rot_eb_e[0] << ","
             << sigma_delta_rot_eb_e[1] << ","
             << sigma_delta_rot_eb_e[2] << "\n";

        return true;
    }

private:
    std::ofstream file;
    bool ok;
};

};

#endif
