#ifndef NAV_EQUATIONS_H
#define NAV_EQUATIONS_H

#include "constants_types.h"
#include "helpers.h"

namespace intnavlib {

    /// \defgroup nav_equations Navigation Equations
    /// Inertial navigation.
    /// @{

    /// \brief Solves the navigation equations in the ECEF frame.
    /// This function propagates the navigation solution (position, velocity, and attitude)
    /// from a previous time step to the current time step using IMU measurements.
    /// The navigation solution update is performed as follows:
    /*!

    \f[
    \mathbf{C}_b^e(+)=\left(\begin{array}{ccc}
    \cos \omega_{i e} \tau_i & \sin \omega_{i e} \tau_i & 0 \\
    -\sin \omega_{i e} \tau_i & \cos \omega_{i e} \tau_i & 0 \\
    0 & 0 & 1
    \end{array}\right) \mathbf{C}_b^e(-) \mathbf{C}_{b+}^{b-} ,
    \f]

    \f[
    \mathbf{f}_{i b}^e=\overline{\mathbf{C}}_b^e \mathbf{f}_{i b}^b, \quad \overline{\mathbf{C}}_b^e=\mathbf{C}_b^e(-) \mathbf{C}_b^{b-}-\frac{1}{2} \boldsymbol{\Omega}_e^e \mathbf{C}_b^e(-) \tau_i ,
    \f]

    \f[
    \mathbf{v}_{e b}^e(+) \approx \mathbf{v}_{e b}^e(-)+\left(\mathbf{f}_{i b}^e+\mathbf{g}_b^e\left(\mathbf{r}_{e b}^e(-)\right)-2 \boldsymbol{\Omega}_{i e}^e \mathbf{v}_{e b}^e(-)\right) \tau_i ,
    \f]

    \f[
    \mathbf{r}_{e b}^e(+)=\mathbf{r}_{e b}^e(-)+\left(\mathbf{v}_{e b}^e(-)+\mathbf{v}_{e b}^e(+)\right) \frac{\tau_i}{2}.
    \f]

    */
    /// \param[in] old_nav The navigation solution at the previous time step.
    /// \param[in] imu_meas IMU measurements (specific force and angular velocity) for the current interval.
    /// \param[in] tor_i The integration time interval (delta_t) in seconds.
    /// \return The updated navigation solution in the ECEF frame.
    NavSolutionEcef navEquationsEcef(const NavSolutionEcef & old_nav, 
                                const ImuMeasurements & imu_meas,
                                const double & tor_i);
    
    /// @}

}

#endif