#ifndef KALMAN_H
#define KALMAN_H

#include <boost/math/distributions/chi_squared.hpp>

#include "constants_types.h"
#include "helpers.h"
#include "nav_equations.h"

namespace intnavlib {

/// \defgroup kalman Kalman
/// Navigation Kalman filter utilities.
/// @{

/// \brief Initializes the error covariance matrix for a navigation Kalman Filter.
/// \param[in] tc_kf_config Configuration parameters for the Kalman Filter,
/// including initial uncertainties for attitude, velocity, position, biases,
/// and clock offset/drift.
/// \return An Eigen::Matrix<double,17,17> representing the initialized error covariance matrix.
Eigen::Matrix<double,17,17> initializePMmatrix(const KfConfig & tc_kf_config);

/// \brief Utility to perform both state prediction and uncertainty propagation (Loose)
StateEstEcef lcPredictKF(const StateEstEcef & state_est_old, 
                            const ImuMeasurements & imu_meas,
                            const KfConfig & kf_config,
                            const double & tor_i);

/// \brief Utility to perform both state prediction and uncertainty propagation (Tight)
StateEstEcef tcPredictKF(const StateEstEcef & state_est_old, 
                            const ImuMeasurements & imu_meas,
                            const KfConfig & kf_config,
                            const double & tor_i);

/// \brief Propagates the error-state uncertainty for a Loosely Coupled (LC) Kalman Filter.
/// This function updates the error covariance matrix based on the system dynamics
/// and noise characteristics over a given integration time.
/*!

Loosely-coupled ECEF Error state:

\f[
\mathbf{x}_{I N S}^e=\left(\begin{array}{c}
\delta \boldsymbol{\Psi}_{e b}^e \\
\delta \mathbf{v}_{e b}^e \\
\delta \mathbf{r}_{e b}^e \\
\mathbf{b}_a \\
\mathbf{b}_g
\end{array}\right)
\f]

Loosely-coupled ECEF Error state first-order linear system model:

\f[
\boldsymbol{\Phi}_{\text {INS }}^e \approx\left[\begin{array}{ccccc}
\mathbf{I}_3-\boldsymbol{\Omega}_e^e \tau_s & 0_3 & 0_3 & 0_3 & \hat{\mathbf{C}}_b^e \tau_s \\
\mathbf{F}_{12}^e \tau_s & \mathbf{I}_3-2 \boldsymbol{\Omega}_e^e \tau_s & \mathbf{F}_{23}^e \tau_s & \hat{\mathbf{C}}_b^e \tau_s & 0_3 \\
0_3 & \mathbf{I}_3 \tau_s & \mathbf{I}_3 & 0_3 & 0_3 \\
0_3 & 0_3 & 0_3 & \mathbf{I}_3 & 0_3 \\
0_3 & 0_3 & 0_3 & 0_3 & \mathbf{I}_3
\end{array}\right]
\f]

\f[
\mathbf{Q}_{I N S}^e \approx  \mathbf{Q'}_{I N S}^e =\left(\begin{array}{ccccc}
S_{r g} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 & 0_3 \\
0_3 & S_{r a} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 \\
0_3 & 0_3 & 0_3 & 0_3 & 0_3 \\
0_3 & 0_3 & 0_3 & S_{b a d} \mathbf{I}_3 & 0_3 \\
0_3 & 0_3 & 0_3 & 0_3 & S_{b g d} \mathbf{I}_3
\end{array}\right) \tau_s
\f]

Loosely-coupled ECEF Error state covariance propagation:

\f[
\mathbf{P}_k^{-} \approx \boldsymbol{\Phi}_{k-1}\left(\mathbf{P}_{k-1}^{+}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime}\right) \boldsymbol{\Phi}_{k-1}^{\mathrm{T}}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime} .
\f]

*/
/// \param[in] P_matrix_old The error covariance matrix from the previous time step.
/// \param[in] old_nav_est_ecef The estimated navigation solution in ECEF frame at the previous time step.
/// \param[in] old_nav_est_ned The estimated navigation solution in NED frame at the previous time step.
/// \param[in] imu_meas IMU measurements used for propagation.
/// \param[in] lc_kf_config Configuration parameters for the LC Kalman Filter,
/// including noise PSDs for gyro, accelerometer, and biases.
/// \param[in] tor_s Integration time (delta_t) in seconds.
/// \return An Eigen::Matrix<double,15,15> representing the propagated error covariance matrix.
Eigen::Matrix<double,15,15> lcPropUnc(const Eigen::Matrix<double,15,15> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & lc_kf_config,
                                        const double & tor_s);

/// \brief Propagates the error-state uncertainty for a Tightly Coupled (TC) Kalman Filter.
/// This function updates the error covariance matrix based on the system dynamics
/// and noise characteristics over a given integration time, including clock states.
/*!

Tightly-coupled ECEF Error state:

\f[
\mathbf{x}_{I N S}^e=\left(\begin{array}{c}
\delta \boldsymbol{\Psi}_{e b}^e \\
\delta \mathbf{v}_{e b}^e \\
\delta \mathbf{r}_{e b}^e \\
\mathbf{b}_a \\
\mathbf{b}_g \\
\delta \rho^a_c \\
\delta \dot{\rho}^a_c
\end{array}\right)
\f]

Tightly-coupled ECEF Error state first-order linear system model:

\f[
\boldsymbol{\Phi}_{\text{INS}}^e \approx \left[\begin{array}{ccccccc}
\mathbf{I}_3 - \boldsymbol{\Omega}_e^e \tau_s & 0_3 & 0_3 & 0_3 & \hat{\mathbf{C}}_b^e \tau_s & 0_{3 \times 1} & 0_{3 \times 1} \\
\mathbf{F}_{12}^e \tau_s & \mathbf{I}_3 - 2\boldsymbol{\Omega}_e^e \tau_s & \mathbf{F}_{23}^e \tau_s & \hat{\mathbf{C}}_b^e \tau_s & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & \mathbf{I}_3 \tau_s & \mathbf{I}_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & 0_3 & 0_3 & \mathbf{I}_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & 0_3 & 0_3 & 0_3 & \mathbf{I}_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 1 & \tau_s \\
0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0 & 1 \\
\end{array}
\right]
\f]

\f[
\mathbf{Q}_{INS}^e \approx \mathbf{Q'}_{I N S}^e =
\left(\begin{array}{ccccccc}
S_{rg} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & S_{ra} \mathbf{I}_3 & 0_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & 0_3 & 0_3 & 0_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & 0_3 & 0_3 & S_{bad} \mathbf{I}_3 & 0_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_3 & 0_3 & 0_3 & 0_3 & S_{bgd} \mathbf{I}_3 & 0_{3 \times 1} & 0_{3 \times 1} \\
0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & S_{cb} & 0 \\
0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0_{1 \times 3} & 0 & S_{cd} \\
\end{array}\right) \tau_s
\f]

Tightly-coupled ECEF Error state covariance propagation:

\f[
\mathbf{P}_k^{-} \approx \boldsymbol{\Phi}_{k-1}\left(\mathbf{P}_{k-1}^{+}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime}\right) \boldsymbol{\Phi}_{k-1}^{\mathrm{T}}+\frac{1}{2} \mathbf{Q}_{k-1}^{\prime} .
\f]
*/
/// \param[in] P_matrix_old The error covariance matrix from the previous time step.
/// \param[in] old_nav_est_ecef The estimated navigation solution in ECEF frame at the previous time step.
/// \param[in] old_nav_est_ned The estimated navigation solution in NED frame at the previous time step.
/// \param[in] imu_meas IMU measurements used for propagation.
/// \param[in] tc_kf_config Configuration parameters for the TC Kalman Filter,
/// including noise PSDs for gyro, accelerometer, biases, and clock states.
/// \param[in] tor_s Integration time (delta_t) in seconds.
/// \return An Eigen::Matrix<double,17,17> representing the propagated error covariance matrix.
Eigen::Matrix<double,17,17> tcPropUnc(const Eigen::Matrix<double,17,17> & P_matrix_old, 
                                        const NavSolutionEcef & old_nav_est_ecef,
                                        const NavSolutionNed & old_nav_est_ned,
                                        const ImuMeasurements & imu_meas,
                                        const KfConfig & tc_kf_config,
                                        const double & tor_s);

/// \brief Performs a Loosely Coupled (LC) Kalman Filter update using an ECEF position measurement.
/// This is a closed-loop error-state KF update that corrects the navigation solution
/// and IMU bias estimates by integrating a position measurement.
/// \param[in] pos_meas The position measurement in ECEF frame.
/// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
/// and its covariance matrix.
/// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
/// \return A StateEstEcefLc structure containing the updated state estimation.
StateEstEcef lcUpdateKFPosEcef (const PosMeasEcef & pos_meas, 
                                    const StateEstEcef & state_est_old,
                                    const double & p_value = 0.99);

/// \brief Performs a Loosely Coupled (LC) Kalman Filter update using a GNSS position and velocity measurement.
/// This is a closed-loop error-state KF update that corrects the navigation solution
/// and IMU bias estimates by integrating a GNSS position and velocity measurement.
/*!

Error state:

\f[
\mathbf{x}_{I N S}^e=\left(\begin{array}{c}
\delta \boldsymbol{\Psi}_{e b}^e \\
\delta \mathbf{v}_{e b}^e \\
\delta \mathbf{r}_{e b}^e \\
\mathbf{b}_a \\
\mathbf{b}_g
\end{array}\right)
\f]

Observation matrix:

\f[
\mathbf{H}_{G, k}^e=\left(\begin{array}{ccccc}
\mathbf{H}_{r 1}^e & 0_3 & -\mathbf{I}_3 & 0_3 & 0_3 \\
\mathbf{H}_{v 1}^e & -\mathbf{I}_3 & 0_3 & 0_3 & \mathbf{H}_{v 5}^e
\end{array}\right)_k
\f]

Innovation:

\f[
\delta \mathbf{z}_{G, k}^{e-} =
\begin{bmatrix}
\hat{\mathbf{r}}_{e a G}^e - \hat{\mathbf{r}}_{e b}^e - \hat{\mathbf{C}}_b^e \mathbf{l}_{b a}^b \\
\hat{\mathbf{v}}_{e a G}^e - \hat{\mathbf{v}}_{e b}^e - \hat{\mathbf{C}}_b^e \left( \hat{\boldsymbol{\omega}}_{i b}^b \wedge \mathbf{l}_{b a}^b \right) + \boldsymbol{\Omega}_e^e \hat{\mathbf{C}}_b^e \mathbf{l}_{b a}^b
\end{bmatrix}_k^{-}
\f]

Error state update:

\f[
\mathbf{K}_k=\mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}\left(\mathbf{H}_k \mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}+\mathbf{R}_k\right)^{-1}
\f]

\f[
\begin{aligned}
\hat{\mathbf{x}}_k^{+} & =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k\left(\mathbf{z}_k-\mathbf{H}_k \hat{\mathbf{x}}_k^{-}\right) \\
& =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k \delta \mathbf{z}_k^{-}
\end{aligned}
\f]

\f[
\mathbf{P}_k^{+}=\left(\mathbf{I}-\mathbf{K}_k \mathbf{H}_k\right) \mathbf{P}_k^{-}
\f]

Navigation solution update

\f[
\begin{gathered}
\hat{\mathbf{C}}_b^\gamma(+)=\delta \hat{\mathbf{C}}_b^{\gamma^{\mathrm{T}}} \hat{\mathbf{C}}_b^\gamma(-) \approx\left(\mathbf{I}_3-\left[\delta \hat{\boldsymbol{\psi}}_{\gamma b}^\gamma \wedge\right]\right) \hat{\mathbf{C}}_b^\gamma(-) \\
\hat{\mathbf{v}}_{\beta b}^\gamma(+)=\hat{\mathbf{v}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{v}}_{\beta b}^\gamma \\
\hat{\mathbf{r}}_{\beta b}^\gamma(+)=\hat{\mathbf{r}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{r}}_{\beta b}^\gamma
\end{gathered}
\f]

*/
/// \param[in] pos_vel_gnss_meas The GNSS position and velocity measurement in ECEF frame.
/// \param[in] state_est_prior The prior estimated state (navigation solution and biases) and its covariance matrix.
/// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
/// \return A StateEstEcefLc structure containing the updated state estimation.
StateEstEcef lcUpdateKFGnssEcef (const GnssPosVelMeasEcef & pos_vel_gnss_meas, 
                                    const StateEstEcef & state_est_old,
                                    const double & p_value = 0.99);

/// \brief Performs a Tightly Coupled (TC) Kalman Filter update using GNSS pseudo-range and pseudo-range rate measurements.
/// This is a closed-loop error-state KF update that corrects the navigation solution,
/// IMU bias estimates, and receiver clock states by integrating raw GNSS measurements.
/*!

Tightly-coupled ECEF Error state:

\f[
\mathbf{x}_{I N S}^e=\left(\begin{array}{c}
\delta \boldsymbol{\Psi}_{e b}^e \\
\delta \mathbf{v}_{e b}^e \\
\delta \mathbf{r}_{e b}^e \\
\mathbf{b}_a \\
\mathbf{b}_g \\
\delta \rho^a_c \\
\delta \dot{\rho}^a_c
\end{array}\right)
\f]

Observation matrix:

\f[
\mathbf{H}_{G, k}^\gamma=\left(\begin{array}{ccccccc}
\frac{\partial \mathbf{z}_\rho}{\partial \delta \boldsymbol{\psi}_{\gamma b}^\gamma} & 0_{m, 3} & \frac{\partial \mathbf{z}_\rho}{\partial \delta \mathbf{r}_{\gamma b}^\gamma} & 0_{m, 3} & 0_{m, 3} & \frac{\partial \mathbf{z}_\rho}{\partial \delta \rho_c^a} & 0_{m, 1} \\
\frac{\partial \mathbf{z}_r}{\partial \delta \boldsymbol{\psi}_{\gamma b}^\gamma} & \frac{\partial \mathbf{z}_r}{\partial \delta \mathbf{v}_{\gamma b}^\gamma} & \frac{\partial \mathbf{z}_r}{\partial \delta \mathbf{r}_{\gamma b}^\gamma} & 0_{m, 3} & \frac{\partial \mathbf{z}_\rho}{\partial \mathbf{b}_g} & 0_{m, 1} & \frac{\partial \mathbf{z}_r}{\partial \delta \dot{\rho}_c^a}
\end{array}\right)_{\mathbf{x}=\hat{\mathbf{x}}_{\bar{k}}}
\f]

Innovation:

\f[
\delta \mathbf{z}_{\bar{G}, k}=\binom{\delta \mathbf{z}_{\rho, k}^{-}}{\delta \mathbf{z}_{r, k}^{-}}, where \quad \delta \mathbf{z}_{\rho, k}^{-}=\left(\tilde{\rho}_{a, C}^1-\hat{\rho}_{a, C}^{1-}, \tilde{\rho}_{a, C}^2-\hat{\rho}_{a, C}^{2-}, \cdots \tilde{\rho}_{a, C}^m-\hat{\rho}_{a, C}^{m-}\right)_k , \quad \delta \mathbf{z}_{r, k}^{-}=\left(\tilde{\dot{\rho}}_{a, C}^1-\hat{\hat{\rho}}_{a, C}^{1-}, \tilde{\hat{\rho}}_{a, C}^2-\hat{\hat{\rho}}_{a, C}^{2-}, \cdots \tilde{\hat{\rho}}_{a, C}^m-\hat{\hat{\rho}}_{a, C}^{m-}\right)_k .
\f]

Error state update:

\f[
\mathbf{K}_k=\mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}\left(\mathbf{H}_k \mathbf{P}_k^{-} \mathbf{H}_k^{\mathrm{T}}+\mathbf{R}_k\right)^{-1}
\f]

\f[
\begin{aligned}
\hat{\mathbf{x}}_k^{+} & =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k\left(\mathbf{z}_k-\mathbf{H}_k \hat{\mathbf{x}}_k^{-}\right) \\
& =\hat{\mathbf{x}}_k^{-}+\mathbf{K}_k \delta \mathbf{z}_k^{-}
\end{aligned}
\f]

\f[
\mathbf{P}_k^{+}=\left(\mathbf{I}-\mathbf{K}_k \mathbf{H}_k\right) \mathbf{P}_k^{-}
\f]

Navigation solution update

\f[
\begin{gathered}
\hat{\mathbf{C}}_b^\gamma(+)=\delta \hat{\mathbf{C}}_b^{\gamma^{\mathrm{T}}} \hat{\mathbf{C}}_b^\gamma(-) \approx\left(\mathbf{I}_3-\left[\delta \hat{\boldsymbol{\psi}}_{\gamma b}^\gamma \wedge\right]\right) \hat{\mathbf{C}}_b^\gamma(-) \\
\hat{\mathbf{v}}_{\beta b}^\gamma(+)=\hat{\mathbf{v}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{v}}_{\beta b}^\gamma \\
\hat{\mathbf{r}}_{\beta b}^\gamma(+)=\hat{\mathbf{r}}_{\beta b}^\gamma(-)-\delta \hat{\mathbf{r}}_{\beta b}^\gamma
\end{gathered}
\f]

*/
/// \param[in] gnss_meas The GNSS measurements, including pseudo-ranges, pseudo-range rates,
/// satellite positions, and satellite velocities.
/// \param[in] state_est_prior The prior estimated state (navigation solution, biases, and clock states)
/// and its covariance matrix.
/// \param[in] tor_s Elapsed time since last update in seconds.
/// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
/// \return A StateEstEcefTc structure containing the updated state estimation.
StateEstEcef tcUpdateKFGnssEcef (const GnssMeasurements & gnss_meas, 
                                    const StateEstEcef & state_est_prior,
                                    const double & tor_s,
                                    const double & p_value = 0.99);

/// \brief Performs a Loosely Coupled (LC) Kalman Filter update using an ECEF position and rotation measurement.
/// This is a closed-loop error-state KF update that corrects the navigation solution
/// and IMU bias estimates by integrating a position and rotation measurement.
/// \param[in] pos_rot_meas The position and rotation measurement in ECEF frame.
/// \param[in] state_est_prior The prior estimated state (navigation solution and biases)
/// and its covariance matrix.
/// \param[in] p_value P-value threshold for Chi squared consistency check: accept update only if residual in the confidence interval.
/// \return A StateEstEcefLc structure containing the updated state estimation.
StateEstEcef lcUpdateKFPosRotEcef (const PosRotMeasEcef & pos_rot_meas, 
                                    const StateEstEcef & state_est_old,
                                    const double & p_value = 0.99);

/// @}

};

#endif