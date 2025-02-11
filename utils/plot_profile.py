import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import contextily as ctx
import geopandas as gpd
from shapely.geometry import Point

# // Column 1: time (sec)
# // Column 2: latitude (deg)
# // Column 3: longitude (deg)
# // Column 4: height (m)
# // Column 5: north velocity (m/s)
# // Column 6: east velocity (m/s)
# // Column 7: down velocity (m/s)
# // Column 8: roll angle of body w.r.t NED (deg)
# // Column 9: pitch angle of body w.r.t NED (deg)
# // Column 10: yaw angle of body w.r.t NED (deg)

def radii_of_curvature(L):
    R_0 = 6378137  # WGS84 Equatorial radius in meters
    e = 0.0818191908425  # WGS84 eccentricity
    temp = 1 - (e * np.sin(L))**2
    R_N = R_0 * (1 - e**2) / temp**1.5
    R_E = R_0 / np.sqrt(temp)
    return R_N, R_E

def plot_profile(csv_file):
    data = pd.read_csv(csv_file)
    time = data.iloc[:, 0].to_numpy()
    latitude = data.iloc[:, 1].to_numpy()
    longitude = data.iloc[:, 2].to_numpy()
    height = data.iloc[:, 3].to_numpy()
    north_velocity = data.iloc[:, 4].to_numpy()
    east_velocity = data.iloc[:, 5].to_numpy()
    down_velocity = data.iloc[:, 6].to_numpy()
    roll_angle = data.iloc[:, 7].to_numpy()
    pitch_angle = data.iloc[:, 8].to_numpy()
    yaw_angle = data.iloc[:, 9].to_numpy()
    
    R_N, R_E = radii_of_curvature(np.deg2rad(latitude.iloc[0]))
    north_displacement = (np.deg2rad(latitude) - np.deg2rad(latitude.iloc[0])) * (R_N + height)
    east_displacement = (np.deg2rad(longitude) - np.deg2rad(longitude.iloc[0])) * (R_N + height) * np.cos(np.radians(latitude.iloc[0]))
    down_displacement = height.iloc[0] - height
    
    fig, axs = plt.subplots(4, 3, figsize=(15, 14))
    fig.tight_layout(pad=5.0)
    
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
    
    # Flattened flight path (latitude vs longitude) with basemap
    geometry = [Point(xy) for xy in zip(longitude, latitude)]
    gdf = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")
    gdf_mercator = gdf.to_crs(epsg=3857)
    
    axs[3, 0].scatter(gdf_mercator.geometry.x, gdf_mercator.geometry.y, color='b', s=5)
    axs[3, 0].set_title('Flattened Flight Path')
    ctx.add_basemap(axs[3, 0], source=ctx.providers.Esri.WorldImagery)
    
    # Convert back to lat/lon for tick labels
    axs[3, 0].set_xticks(gdf_mercator.geome
try.y[::len(gdf_mercator)//5])
    axs[3, 0].set_yticklabels([f'{lat:.2f}' for lat in latitude[::len(latitude)//5]])
    axs[3, 0].set_xlabel('Longitude (deg)')
    axs[3, 0].set_ylabel('Latitude (deg)')
    
    # Altitude profile (height vs time)
    axs[3, 1].plot(time, height, color='g', linewidth=1.5)
    axs[3, 1].set_title('Altitude Profile')
    axs[3, 1].set_xlabel('Time (s)')
    axs[3, 1].set_ylabel('Height (m)')
    axs[3, 1].grid()
    
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
