#ifndef MOTION_PROFILE_WRITER_H
#define MOTION_PROFILE_WRITER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "helpers.h"

using namespace helpers;

// Column 1: time (sec)
// Column 2: latitude (deg)
// Column 3: longitude (deg)
// Column 4: height (m)
// Column 5: north velocity (m/s)
// Column 6: east velocity (m/s)
// Column 7: down velocity (m/s)
// Column 8: roll angle of body w.r.t NED (deg)
// Column 9: pitch angle of body w.r.t NED (deg)
// Column 10: yaw angle of body w.r.t NED (deg)

class MotionProfileWriter {
public:
    MotionProfileWriter(const std::string& filename) : file(filename), ok(false) {
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Write header
        // file << "time,latitude,longitude,height,vx,vy,vz,roll,pitch,yaw\n";
        ok = true;
    }


    bool writeNextRow(const NavSolutionNed& row) {
        if (!file.good()) return false;

        // Convert latitude and longitude from radians to degrees
        double latitude = helpers::radToDeg(row.latitude);
        double longitude = helpers::radToDeg(row.longitude);

        // Extract velocity components
        Eigen::Vector3d v_eb_n = row.v_eb_n;

        // Convert rotation matrix to roll, pitch, yaw (in radians)
        Eigen::Vector3d rpy = helpers::rToRpy(row.C_b_n.transpose());

        // Convert roll, pitch, yaw to degrees
        double roll = helpers::radToDeg(rpy[0]);
        double pitch = helpers::radToDeg(rpy[1]);
        double yaw = helpers::radToDeg(rpy[2]);

        // Write the row to the CSV file
        file << row.time << ","
             << latitude << ","
             << longitude << ","
             << row.height << ","
             << v_eb_n[0] << ","
             << v_eb_n[1] << ","
             << v_eb_n[2] << ","
             << roll << ","
             << pitch << ","
             << yaw << "\n";

        return true;
    }

private:
    std::ofstream file;
    bool ok;
};

#endif
