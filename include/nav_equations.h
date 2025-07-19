#ifndef NAV_EQUATIONS_H
#define NAV_EQUATIONS_H

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    /// \brief Solves the navigation equations in the ECEF frame.
    /// This function propagates the navigation solution (position, velocity, and attitude)
    /// from a previous time step to the current time step using IMU measurements.
    /// \param[in] old_nav The navigation solution at the previous time step.
    /// \param[in] imu_meas IMU measurements (specific force and angular velocity) for the current interval.
    /// \param[in] tor_i The integration time interval (delta_t) in seconds.
    /// \return The updated navigation solution in the ECEF frame.
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i);

}

#endif