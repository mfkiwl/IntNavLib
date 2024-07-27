import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_errors(csv_file):
    # Read the CSV file into a DataFrame
    data = pd.read_csv(csv_file)
    
    # Extract the relevant columns
    time = data.iloc[:, 0]

    x_pos_error = abs(data.iloc[:, 1])
    y_pos_error = abs(data.iloc[:, 2])
    z_pos_error = abs(data.iloc[:, 3])
    x_vel_error = abs(data.iloc[:, 4])
    y_vel_error = abs(data.iloc[:, 5])
    z_vel_error = abs(data.iloc[:, 6])
    roll_error = abs(data.iloc[:, 7])
    pitch_error = abs(data.iloc[:, 8])
    yaw_error = abs(data.iloc[:, 9])

    x_pos_error_sigma = data.iloc[:, 10] * 3
    y_pos_error_sigma = data.iloc[:, 11] * 3
    z_pos_error_sigma = data.iloc[:, 12] * 3
    x_vel_error_sigma = data.iloc[:, 13] * 3
    y_vel_error_sigma = data.iloc[:, 14] * 3
    z_vel_error_sigma = data.iloc[:, 15] * 3
    roll_error_sigma = data.iloc[:, 16] * 3
    pitch_error_sigma = data.iloc[:, 17] * 3
    yaw_error_sigma = data.iloc[:, 18] * 3

    # Create a 3x3 grid of subplots
    fig, axs = plt.subplots(3, 3, figsize=(15, 10))
    fig.tight_layout(pad=5.0)
    
    # Plot each error on the appropriate subplot
    axs[0, 0].plot(time, x_pos_error, color=[0.9, 0.45, 0], linewidth=1.5)
    axs[0, 0].plot(time, x_pos_error_sigma, color=[0.9, 0.45, 0], linestyle="dashed", linewidth=1.5)
    axs[0, 0].set_title('X position error, m')
    
    axs[0, 1].plot(time, y_pos_error, color=[0, 0.9, 0.45], linewidth=1.5)
    axs[0, 1].plot(time, y_pos_error_sigma, color=[0, 0.9, 0.45], linestyle="dashed", linewidth=1.5)
    axs[0, 1].set_title('Y position error, m')
    
    axs[0, 2].plot(time, z_pos_error, color=[0.45, 0, 0.9], linewidth=1.5)
    axs[0, 2].plot(time, z_pos_error_sigma, color=[0.45, 0, 0.9], linestyle="dashed", linewidth=1.5)
    axs[0, 2].set_title('Z position error, m')
    
    axs[1, 0].plot(time, x_vel_error, color=[0.9, 0, 0.45], linewidth=1.5)
    axs[1, 0].plot(time, x_vel_error_sigma, color=[0.9, 0, 0.45], linestyle="dashed", linewidth=1.5)
    axs[1, 0].set_title('X velocity error, m/s')
    
    axs[1, 1].plot(time, y_vel_error, color=[0.45, 0.9, 0], linewidth=1.5)
    axs[1, 1].plot(time, y_vel_error_sigma, color=[0.45, 0.9, 0], linestyle="dashed", linewidth=1.5)
    axs[1, 1].set_title('Y velocity error, m/s')
    
    axs[1, 2].plot(time, z_vel_error, color=[0, 0.45, 0.9], linewidth=1.5)
    axs[1, 2].plot(time, z_vel_error_sigma, color=[0, 0.45, 0.9], linestyle="dashed", linewidth=1.5)
    axs[1, 2].set_title('Z velocity error, m/s')
    
    axs[2, 0].plot(time, roll_error_sigma, color=[0, 0.7, 0.7], linestyle="dashed", linewidth=1.5)
    axs[2, 0].plot(time, roll_error, color=[0, 0.7, 0.7], linewidth=1.5)
    axs[2, 0].set_title('Roll error, deg')
    axs[2, 0].set_xlabel('Time, s')
    
    axs[2, 1].plot(time, pitch_error, color=[0.7, 0, 0.7], linewidth=1.5)
    axs[2, 1].plot(time, pitch_error_sigma, color=[0.7, 0, 0.7], linestyle="dashed", linewidth=1.5)
    axs[2, 1].set_title('Pitch error, deg')
    axs[2, 1].set_xlabel('Time, s')
    
    axs[2, 2].plot(time, yaw_error, color=[0.7, 0.7, 0], linewidth=1.5)
    axs[2, 2].plot(time, yaw_error_sigma, color=[0.7, 0.7, 0], linestyle="dashed", linewidth=1.5)
    axs[2, 2].set_title('Yaw error, deg')
    axs[2, 2].set_xlabel('Time, s')
    
    # Display the plot
    plt.show()


def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_errors.py <csv_file>")
        return
    
    csv_file = sys.argv[1]
    plot_errors(csv_file)

if __name__ == "__main__":
    main()
