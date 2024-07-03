#ifndef MOTION_PROFILE_READER_H
#define MOTION_PROFILE_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "helpers.h"

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

using namespace helpers;

class MotionProfileReader {
public:
    MotionProfileReader(const std::string& filename) : file(filename), ok(false) {
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        // Verify the column count by reading the first row
        std::string header;
        if (std::getline(file, header)) {
            std::stringstream ss(header);
            std::string token;
            int column_count = 0;
            while (std::getline(ss, token, ',')) {
                column_count++;
            }

            if (column_count == 10) { // Ensure 10 columns as expected
                ok = true;
            } else {
                std::cerr << "Input file has the wrong number of columns: " << column_count << std::endl;
            }
        }
    }

    bool isOk() const {
        return ok;
    }

    bool readNextRow(NavSolutionNed& row) {
        if (!file.good()) return false;

        std::string line;
        if (std::getline(file, line)) {
            std::stringstream ss(line);
            std::vector<double> values;
            std::string token;

            while (std::getline(ss, token, ',')) {
                values.push_back(std::stod(token));
            }

            if (values.size() == 10) {
                row.time = values[0];
                row.latitude = helpers::degToRad(values[1]);
                row.longitude = helpers::degToRad(values[2]);
                row.height = values[3];
                row.v_eb_n = Eigen::Vector3d(values[4], values[5], values[6]);

                double roll = helpers::degToRad(values[7]);
                double pitch = helpers::degToRad(values[8]);
                double yaw = helpers::degToRad(values[9]);

                Eigen::Vector3d rpy;
                rpy << roll, pitch, yaw;
                row.C_b_n = helpers::rpyToR(rpy);

                return true;
            }
        }
        return false;
    }

private:
    std::ifstream file;
    bool ok;
};


#endif