#include <thread>
#include <chrono>

#include "helpers.h"
#include "lc_ins_gnss.h"
#include "motion_profile_reader.h"

int main(int argc, char** argv)
{   

    // Init data loader (reads motion profile, generates measurements)

    // Init filter

    // Iterate data loader

    MotionProfileReader reader("profile.csv");

    MotionProfileRow row;
    while (reader.readNextRow(row)) {
        std::cout << "Time: " << row.time << ", Latitude (radians): " << row.latitude
                  << ", Longitude (radians): " << row.longitude << std::endl;
        // Process other fields as necessary
    }

    return 0;
}