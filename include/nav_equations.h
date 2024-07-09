#ifndef NAV_EQUATIONS_H
#define NAV_EQUATIONS_H

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    // Nav equations in ECEF
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i);

}

#endif