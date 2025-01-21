# Custom GnssPosVelMeasEcef.msg to match my struct

std_msgs/Header header

# Position vector in ECEF (Earth-Centered, Earth-Fixed) coordinates
geometry_msgs/Vector3 r_ea_e

# Velocity vector in ECEF coordinates
geometry_msgs/Vector3 v_ea_e

# Covariance matrix (flattened to a 1D array for compatibility)
# float64[36] cov_mat