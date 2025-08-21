import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_errors(csv_file):
    # Read the CSV file into a DataFrame
    data = pd.read_csv(csv_file)
    
    # Extract the relevant columns
    time = data.iloc[:, 0].to_numpy()

    x_pos_error = data.iloc[:, 1].to_numpy()
    y_pos_error = data.iloc[:, 2].to_numpy()
    z_pos_error = data.iloc[:, 3].to_numpy()
    x_vel_error = data.iloc[:, 4].to_numpy()
    y_vel_error = data.iloc[:, 5].to_numpy()
    z_vel_error = data.iloc[:, 6].to_numpy()
    roll_error = data.iloc[:, 7].to_numpy()
    pitch_error = data.iloc[:, 8].to_numpy()
    yaw_error = data.iloc[:, 9].to_numpy()

    # sigmas (±3σ)
    x_pos_error_sigma = data.iloc[:, 10].to_numpy() * 3
    y_pos_error_sigma = data.iloc[:, 11].to_numpy() * 3
    z_pos_error_sigma = data.iloc[:, 12].to_numpy() * 3
    x_vel_error_sigma = data.iloc[:, 13].to_numpy() * 3
    y_vel_error_sigma = data.iloc[:, 14].to_numpy() * 3
    z_vel_error_sigma = data.iloc[:, 15].to_numpy() * 3
    roll_error_sigma = data.iloc[:, 16].to_numpy() * 3
    pitch_error_sigma = data.iloc[:, 17].to_numpy() * 3
    yaw_error_sigma = data.iloc[:, 18].to_numpy() * 3

    # Create a 3x3 grid of subplots
    fig, axs = plt.subplots(3, 3, figsize=(15, 10))
    fig.tight_layout(pad=5.0)
    
    # Utility function to plot error with ±sigma bounds
    def plot_with_bounds(ax, time, error, sigma, color, title, xlabel=None):
        ax.plot(time, error, color=color, linewidth=1.5)
        ax.plot(time, sigma, color=color, linestyle="dashed", linewidth=1.0)
        ax.plot(time, -sigma, color=color, linestyle="dashed", linewidth=1.0)
        ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
    
    # Position errors
    plot_with_bounds(axs[0, 0], time, x_pos_error, x_pos_error_sigma, [0.9, 0.45, 0], 'X position error, m')
    plot_with_bounds(axs[0, 1], time, y_pos_error, y_pos_error_sigma, [0, 0.9, 0.45], 'Y position error, m')
    plot_with_bounds(axs[0, 2], time, z_pos_error, z_pos_error_sigma, [0.45, 0, 0.9], 'Z position error, m')
    
    # Velocity errors
    plot_with_bounds(axs[1, 0], time, x_vel_error, x_vel_error_sigma, [0.9, 0, 0.45], 'X velocity error, m/s')
    plot_with_bounds(axs[1, 1], time, y_vel_error, y_vel_error_sigma, [0.45, 0.9, 0], 'Y velocity error, m/s')
    plot_with_bounds(axs[1, 2], time, z_vel_error, z_vel_error_sigma, [0, 0.45, 0.9], 'Z velocity error, m/s')
    
    # Rotation errors
    plot_with_bounds(axs[2, 0], time, roll_error, roll_error_sigma, [0, 0.7, 0.7], 'X rotation error, deg', xlabel='Time, s')
    plot_with_bounds(axs[2, 1], time, pitch_error, pitch_error_sigma, [0.7, 0, 0.7], 'Y rotation error, deg', xlabel='Time, s')
    plot_with_bounds(axs[2, 2], time, yaw_error, yaw_error_sigma, [0.7, 0.7, 0], 'Z rotation error, deg', xlabel='Time, s')
    
    # Display the plot
    plt.show()


def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_errors_sigmas.py <csv_file>")
        return
    
    csv_file = sys.argv[1]
    plot_errors(csv_file)

if __name__ == "__main__":
    main()
