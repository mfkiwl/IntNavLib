#ifndef HELPERS_H
#define HELPERS_H

#include <random>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>

#include "constants_types.h"

namespace intnavlib {

    /// \defgroup Helpers
    /// @{
    
    /// \brief Gets gravity vector at given ECEF position, in ECEF frame.
    /// \param[in] r_eb_e Position vector in ECEF frame.
    /// \return Gravity vector in ECEF frame.
    Eigen::Vector3d gravityEcef(const Eigen::Vector3d & r_eb_e);

    /// \brief Converts a navigation solution from NED (North, East, Down) frame to ECEF (Earth-Centered, Earth-Fixed) frame.
    /// \param[in] nav_sol_ned Navigation solution in NED frame.
    /// \return Navigation solution in ECEF frame.
    NavSolutionEcef nedToEcef(const NavSolutionNed & nav_sol_ned);

    /// \brief Converts a navigation solution from ECEF (Earth-Centered, Earth-Fixed) frame to NED (North, East, Down) frame.
    /// \param[in] nav_sol_ecef Navigation solution in ECEF frame.
    /// \return Navigation solution in NED frame.
    NavSolutionNed ecefToNed(const NavSolutionEcef & nav_sol_ecef);

    /// \brief Compute errors from state estimate and ground truth.
    /// \param[in] state_est_ecef State estimate.
    /// \param[in] true_nav_ecef grounud truth.
    /// \return Errors and standard deviations.
    ErrorsSigmasEcef getErrorsSigmasEcef(const StateEstEcef & state_est_ecef, const NavSolutionEcef & true_nav_ecef);

    /// \brief Converts degrees to radians.
    /// \param[in] degrees Angle in degrees.
    /// \return Angle in radians.
    inline double degToRad(const double & degrees) {
        return degrees * 0.01745329252;
    }

    /// \brief Converts radians to degrees.
    /// \param[in] rads Angle in radians.
    /// \return Angle in degrees.
    inline double radToDeg(const double & rads) {
        return rads / 0.01745329252;
    }

    /// \brief Converts Euler angles (roll, pitch, yaw) to a rotation matrix C.
    /// The Euler sequence is the following: C = Rz(yaw) * Ry(pitch) * Rx(roll).
    /// \param[in] rpy Eigen::Vector3d containing roll, pitch, and yaw angles (in radians).
    /// \return 3x3 rotation matrix.
    Eigen::Matrix3d eulerToDcm(const Eigen::Vector3d & rpy);

    /// \brief Converts a rotation matrix to Euler angles (roll, pitch, yaw).
    /// The Euler sequence is the following: R = Rz(yaw) * Ry(pitch) * Rx(roll).
    /// \param[in] C 3x3 rotation matrix.
    /// \return Eigen::Vector3d containing roll, pitch, and yaw angles (in radians).
    Eigen::Vector3d dcmToEuler(const Eigen::Matrix3d & C);

    /// \brief Gets the skew-symmetric matrix from a 3D vector.
    /// For a vector `a = [ax, ay, az]`, the skew-symmetric matrix `S` is:
    /// `[ 0  -az  ay ]`
    /// `[ az  0  -ax ]`
    /// `[-ay  ax  0  ]`
    /// \param[in] a Input 3D vector.
    /// \return 3x3 skew-symmetric matrix.
    Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d & a);

    /// \brief Signum function.
    /// Returns -1 if `val < 0`, 0 if `val == 0`, and 1 if `val > 0`.
    /// \tparam T Type of the input value.
    /// \param[in] val Input value.
    /// \return Sign of the input value.
    template <typename T> inline T sgn(const T & val) {
    return T((T(0) < val) - (val < T(0)));
    }

    /// \brief Gets the current date and time as a formatted string.
    /// The format is "YYYY-MM-DD_HH-MM-SS".
    /// \return A string representing the current date and time.
    std::string getCurrentDateTime();

    /// \brief Gets the radii of curvature (meridian and transverse) from latitude.
    /// \param[in] L Latitude in radians.
    /// \return Eigen::Vector2d where x is the meridian radius and y is the transverse radius.
    Eigen::Vector2d radiiOfCurvature(double L);

    /// \brief Calculates position, velocity, and attitude errors in the local NED (North, East, Down) frame.
    /// \param[in] true_nav_sol True navigation solution in NED frame.
    /// \param[in] est_nav_sol Estimated navigation solution in NED frame.
    /// \return Structure containing the calculated errors.
    ErrorsNed calculateErrorsNed(const NavSolutionNed & true_nav_sol, 
                                const NavSolutionNed & est_nav_sol);

    /// \brief Extracts the vector from a skew-symmetric matrix.
    /// This is the inverse operation of `skewSymmetric`.
    /// For a skew-symmetric matrix `S`:
    /// `[ 0  -az  ay ]`
    /// `[ az  0  -ax ]`
    /// `[-ay  ax  0  ]`
    /// the function returns the vector `[ax, ay, az]`.
    /// \param[in] S 3x3 skew-symmetric matrix.
    /// \return The 3D vector from which the skew-symmetric matrix was formed.
    Eigen::Vector3d deSkew(const Eigen::Matrix3d & S);
    
    /// \brief Get IMU errors for a tactical grade IMU.
    ImuErrors tacticalImuErrors();
    /// \brief Get default GNSS config.
    GnssConfig defaultGnssConfig();
    /// \brief Get default KfConfig.
    KfConfig tacticalImuKFConfig(); 
    
    /// \brief Text file writer utility
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
            Eigen::Vector3d delta_rot_eb_e = row.delta_rot_eb_e * kRadToDeg;
            Eigen::Vector3d sigma_delta_rot_eb_e = row.sigma_delta_rot_eb_e * kRadToDeg;
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

    /// \brief Profile reader utility
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
    class MotionProfileReader {
        public:
            MotionProfileReader(const std::string& filename) : file(filename), ok(false) {
                if (!file.is_open()) {
                    std::cerr << "Failed to open file: " << filename << std::endl;
                    return;
                }
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
                        row.latitude = degToRad(values[1]);
                        row.longitude = degToRad(values[2]);
                        row.height = values[3];
                        row.v_eb_n = Eigen::Vector3d(values[4], values[5], values[6]);
                        double roll = degToRad(values[7]);
                        double pitch = degToRad(values[8]);
                        double yaw = degToRad(values[9]);
                        Eigen::Vector3d rpy;
                        rpy << roll, pitch, yaw;
                        row.C_b_n = eulerToDcm(rpy);

                        return true;
                    }
                    else{
                        throw std::runtime_error("Wrong number of columns in input profile.");
                    }
                }
                return false;
            }
        private:
            std::ifstream file;
            bool ok;
    };

    /// \brief Profile writer utility
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
                // Since we will use this csv file to compute errors between profiles, its important to set precision
                file << std::fixed << std::setprecision(20);
                // Write header
                // file << "time,latitude,longitude,height,vx,vy,vz,roll,pitch,yaw\n";
                ok = true;
            }
            bool writeNextRow(const NavSolutionNed& row) {
                if (!file.good()) return false;
                // Convert latitude and longitude from radians to degrees
                double latitude = radToDeg(row.latitude);
                double longitude = radToDeg(row.longitude);
                double height = row.height;
                // Extract velocity components
                Eigen::Vector3d v_eb_n = row.v_eb_n;
                // Convert rotation matrix to roll, pitch, yaw (in radians)
                Eigen::Vector3d rpy = dcmToEuler(row.C_b_n);
                // Convert roll, pitch, yaw to degrees
                double roll = radToDeg(rpy[0]);
                double pitch = radToDeg(rpy[1]);
                double yaw = radToDeg(rpy[2]);
                // Write the row to the CSV file
                file << row.time << ","
                    << latitude << ","
                    << longitude << ","
                    << height << ","
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

    /// \brief NED Errors writer
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
                Eigen::Vector3d delta_rot_nb_n = row.delta_rot_nb_n * kRadToDeg;

                // Write the row to the CSV file
                file << row.time << ","
                    << delta_r_eb_n[0] << "," // north position error
                    << delta_r_eb_n[1] << "," // east position error
                    << delta_r_eb_n[2] << "," // down position error
                    << delta_v_eb_n[0] << "," // north velocity error
                    << delta_v_eb_n[1] << "," // east velocity error
                    << delta_v_eb_n[2] << "," // down velocity error
                    << delta_rot_nb_n[0] << "," // roll error
                    << delta_rot_nb_n[1] << "," // pitch error
                    << delta_rot_nb_n[2] << "\n"; // yaw error

                return true;
            }

        private:
            std::ofstream file;
            bool ok;
        };
                                
    /// @}                      
                                
};

#endif