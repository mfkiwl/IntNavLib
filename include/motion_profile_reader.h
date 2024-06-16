#ifndef MOTION_PROFILE_READER_H
#define MOTION_PROFILE_READER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "helpers.h"

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

            if (column_count == 10) {
                ok = true;
            } else {
                std::cerr << "Input file has the wrong number of columns" << std::endl;
            }
        }
    }

    bool isOk() const {
        return ok;
    }

    bool readNextRow(MotionProfileRow& row) {
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
                row.latitude = degToRad(values[1]);
                row.longitude = degToRad(values[2]);
                row.height = values[3];
                row.north_velocity = values[4];
                row.east_velocity = values[5];
                row.down_velocity = values[6];
                row.roll_angle = degToRad(values[7]);
                row.pitch_angle = degToRad(values[8]);
                row.yaw_angle = degToRad(values[9]);
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