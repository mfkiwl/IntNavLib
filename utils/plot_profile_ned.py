import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import contextily as ctx
import geopandas as gpd
from shapely.geometry import Point

def radii_of_curvature(L):
    kR0 = 6378137.0  # WGS84 Equatorial radius [m]
    e = 0.0818191908425  # eccentricity
    temp = 1 - (e * np.sin(L)) ** 2
    R_N = kR0 * (1 - e ** 2) / temp ** 1.5
    R_E = kR0 / np.sqrt(temp)
    return R_N, R_E

def plot_profile(csv_file):
    data = pd.read_csv(csv_file, header=None)
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

    # Compute displacements
    R_N, R_E = radii_of_curvature(np.deg2rad(latitude))
    north_displacement = (np.deg2rad(latitude) - np.deg2rad(latitude[0])) * (R_N + height)
    east_displacement = (np.deg2rad(longitude) - np.deg2rad(longitude[0])) * (R_E + height) * np.cos(np.deg2rad(latitude))
    down_displacement = height - height[0]

    fig, axs = plt.subplots(4, 3, figsize=(15, 14))
    fig.tight_layout(pad=5.0)

    axs[0, 0].plot(time, north_displacement, color=[0.9, 0.45, 0], lw=1.5)
    axs[0, 0].set_title("North displacement [m]")

    axs[0, 1].plot(time, east_displacement, color=[0, 0.9, 0.45], lw=1.5)
    axs[0, 1].set_title("East displacement [m]")

    axs[0, 2].plot(time, down_displacement, color=[0.45, 0, 0.9], lw=1.5)
    axs[0, 2].set_title("Down displacement [m]")

    axs[1, 0].plot(time, north_velocity, color=[0.9, 0, 0.45], lw=1.5)
    axs[1, 0].set_title("North velocity [m/s]")

    axs[1, 1].plot(time, east_velocity, color=[0.45, 0.9, 0], lw=1.5)
    axs[1, 1].set_title("East velocity [m/s]")

    axs[1, 2].plot(time, down_velocity, color=[0, 0.45, 0.9], lw=1.5)
    axs[1, 2].set_title("Down velocity [m/s]")

    axs[2, 0].plot(time, roll_angle, color=[0, 0.7, 0.7], lw=1.5)
    axs[2, 0].set_title("Bank [deg]")
    axs[2, 0].set_xlabel("Time [s]")

    axs[2, 1].plot(time, pitch_angle, color=[0.7, 0, 0.7], lw=1.5)
    axs[2, 1].set_title("Elevation [deg]")
    axs[2, 1].set_xlabel("Time [s]")

    axs[2, 2].plot(time, yaw_angle, color=[0.7, 0.7, 0], lw=1.5)
    axs[2, 2].set_title("Heading [deg]")
    axs[2, 2].set_xlabel("Time [s]")

    # Flattened flight path with basemap
    geometry = [Point(xy) for xy in zip(longitude, latitude)]
    gdf = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326").to_crs(epsg=3857)

    axs[3, 0].scatter(gdf.geometry.x, gdf.geometry.y, c="b", s=5)
    axs[3, 0].set_title("Flattened Flight Path")
    ctx.add_basemap(axs[3, 0], source=ctx.providers.Esri.WorldImagery)

    # Label approximate geographic ticks
    sample_idx = np.linspace(0, len(latitude) - 1, 5, dtype=int)
    xticks = gdf.geometry.x.iloc[sample_idx]
    yticks = gdf.geometry.y.iloc[sample_idx]
    axs[3, 0].set_xticks(xticks)
    axs[3, 0].set_xticklabels([f"{longitude[i]:.2f}" for i in sample_idx])
    axs[3, 0].set_yticks(yticks)
    axs[3, 0].set_yticklabels([f"{latitude[i]:.2f}" for i in sample_idx])
    axs[3, 0].set_xlabel("Longitude [deg]")
    axs[3, 0].set_ylabel("Latitude [deg]")

    # Altitude vs time
    axs[3, 1].plot(time, height, color="g", lw=1.5)
    axs[3, 1].set_title("Altitude Profile")
    axs[3, 1].set_xlabel("Time [s]")
    axs[3, 1].set_ylabel("Height [m]")
    axs[3, 1].grid(True)

    # Hide unused subplot (bottom right)
    axs[3, 2].axis("off")

    plt.show()

def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_profile_ned.py <csv_file>")
        return
    plot_profile(sys.argv[1])

if __name__ == "__main__":
    main()
