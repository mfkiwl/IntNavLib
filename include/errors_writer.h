#ifndef ERRORS_WRITER_H
#define ERRORS_WRITER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "helpers.h"

using namespace intnavlib;

// Column 1: time (sec)
// Column 2: north position error (m)
// Column 3: east position error (m)
// Column 4: down position error (m)
// Column 5: north velocity error (m/s)
// Column 6: east velocity error (m/s)
// Column 7: down velocity error (m/s)
// Column 8: roll component of NED attitude error (deg)
// Column 9: pitch component of NED attitude error (deg)
// Column 10: yaw component of NED attitude error (deg)

class ErrorsWriter {
public:
    ErrorsWriter(const std::string& filename) : file(filename), ok(false) {
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Write header
        // file << "time,d_north,d_east,d_down,d_v_north,d_v_east,d_v_down,d_roll,d_pitch,d_yaw\n";
        ok = true;
    }

    bool writeNextRow(const ErrorsNed& row) {
        if (!file.good()) return false;

         // Extract position errors (in meters)
        Eigen::Vector3d delta_r_eb_n = row.delta_r_eb_n;

        // Extract velocity errors (in m/s)
        Eigen::Vector3d delta_v_eb_n = row.delta_v_eb_n;

        // Convert Euler angles (in radians) to degrees
        Eigen::Vector3d delta_eul_nb_n = row.delta_eul_nb_n * rad_to_deg;

        // Write the row to the CSV file
        file << row.time << ","
             << delta_r_eb_n[0] << "," // north position error
             << delta_r_eb_n[1] << "," // east position error
             << delta_r_eb_n[2] << "," // down position error
             << delta_v_eb_n[0] << "," // north velocity error
             << delta_v_eb_n[1] << "," // east velocity error
             << delta_v_eb_n[2] << "," // down velocity error
             << delta_eul_nb_n[0] << "," // roll error
             << delta_eul_nb_n[1] << "," // pitch error
             << delta_eul_nb_n[2] << "\n"; // yaw error

        return true;
    }

private:
    std::ofstream file;
    bool ok;
};

#endif
