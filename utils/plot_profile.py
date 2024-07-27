import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Column 1: time (sec)
# Column 2: latitude (deg)
# Column 3: longitude (deg)
# Column 4: height (m)
# Column 5: north velocity (m/s)
# Column 6: east velocity (m/s)
# Column 7: down velocity (m/s)
# Column 8: roll angle of body w.r.t NED (deg)
# Column 9: pitch angle of body w.r.t NED (deg)
# Column 10: yaw angle of body w.r.t NED (deg)

def radii_of_curvature(L):

    # Parameters
    R_0 = 6378137  # WGS84 Equatorial radius in meters
    e = 0.0818191908425  # WGS84 eccentricity
    
    # Calculate meridian radius of curvature
    temp = 1 - (e * np.sin(L))**2
    R_N = R_0 * (1 - e**2) / temp**1.5
    
    # Calculate transverse radius of curvature
    R_E = R_0 / np.sqrt(temp)
    
    return R_N, R_E

def plot_profile(csv_file):
    # Read the CSV file into a DataFrame
    data = pd.read_csv(csv_file)
    
    # Extract the relevant columns
    time = data.iloc[:, 0]
    latitude = np.deg2rad(data.iloc[:, 1])
    longitude = np.deg2rad(data.iloc[:, 2])
    height = data.iloc[:, 3]
    north_velocity = data.iloc[:, 4]
    east_velocity = data.iloc[:, 5]
    down_velocity = data.iloc[:, 6]
    roll_angle = data.iloc[:, 7]
    pitch_angle = data.iloc[:, 8]
    yaw_angle = data.iloc[:, 9]
    
    # Calculate radii of curvature
    R_N, R_E = radii_of_curvature(latitude.iloc[0])
    
    # Calculate displacements
    north_displacement = (latitude - latitude.iloc[0]) * (R_N + height)
    east_displacement = (longitude - longitude.iloc[0]) * (R_N + height) * np.cos(np.radians(latitude.iloc[0]))
    down_displacement = height.iloc[0] - height
    
    # Create a 3x3 grid of subplots
    fig, axs = plt.subplots(3, 3, figsize=(15, 10))
    fig.tight_layout(pad=5.0)
    
    # Plot each profile on the appropriate subplot
    axs[0, 0].plot(time, north_displacement, color=[0.9, 0.45, 0], linewidth=1.5)
    axs[0, 0].set_title('North displacement, m')
    
    axs[0, 1].plot(time, east_displacement, color=[0, 0.9, 0.45], linewidth=1.5)
    axs[0, 1].set_title('East displacement, m')
    
    axs[0, 2].plot(time, down_displacement, color=[0.45, 0, 0.9], linewidth=1.5)
    axs[0, 2].set_title('Down displacement, m')
    
    axs[1, 0].plot(time, north_velocity, color=[0.9, 0, 0.45], linewidth=1.5)
    axs[1, 0].set_title('North velocity, m/s')
    
    axs[1, 1].plot(time, east_velocity, color=[0.45, 0.9, 0], linewidth=1.5)
    axs[1, 1].set_title('East velocity, m/s')
    
    axs[1, 2].plot(time, down_velocity, color=[0, 0.45, 0.9], linewidth=1.5)
    axs[1, 2].set_title('Down velocity, m/s')
    
    axs[2, 0].plot(time, roll_angle, color=[0, 0.7, 0.7], linewidth=1.5)
    axs[2, 0].set_title('Bank, deg')
    axs[2, 0].set_xlabel('Time, s')
    
    axs[2, 1].plot(time, pitch_angle, color=[0.7, 0, 0.7], linewidth=1.5)
    axs[2, 1].set_title('Elevation, deg')
    axs[2, 1].set_xlabel('Time, s')
    
    axs[2, 2].plot(time, yaw_angle, color=[0.7, 0.7, 0], linewidth=1.5)
    axs[2, 2].set_title('Heading, deg')
    axs[2, 2].set_xlabel('Time, s')
    
    # Display the plot
    plt.show()

def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_profile.py <csv_file>")
        return
    
    csv_file = sys.argv[1]
    plot_profile(csv_file)

if __name__ == "__main__":
    main()
