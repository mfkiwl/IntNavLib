import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Column 1: time (sec)
# Column 2: north position error (m)
# Column 3: east position error (m)
# Column 4: down position error (m)
# Column 5: north velocity error (m/s)
# Column 6: east velocity error (m/s)
# Column 7: down velocity error (m/s)
# Column 8: roll component of NED attitude error (deg)
# Column 9: pitch component of NED attitude error (deg)
# Column 10: yaw component of NED attitude error (deg)

def plot_errors(csv_file):
    # Read the CSV file into a DataFrame
    data = pd.read_csv(csv_file)
    
    # Extract the relevant columns
    time = data.iloc[:, 0]
    north_pos_error = data.iloc[:, 1]
    east_pos_error = data.iloc[:, 2]
    down_pos_error = data.iloc[:, 3]
    north_vel_error = data.iloc[:, 4]
    east_vel_error = data.iloc[:, 5]
    down_vel_error = data.iloc[:, 6]
    roll_error = data.iloc[:, 7]
    pitch_error = data.iloc[:, 8]
    yaw_error = data.iloc[:, 9]
    
    # Create a 3x3 grid of subplots
    fig, axs = plt.subplots(3, 3, figsize=(15, 10))
    fig.tight_layout(pad=5.0)
    
    # Plot each error on the appropriate subplot
    axs[0, 0].plot(time, north_pos_error, color=[0.9, 0.45, 0], linewidth=1.5)
    axs[0, 0].set_title('North position error, m')
    
    axs[0, 1].plot(time, east_pos_error, color=[0, 0.9, 0.45], linewidth=1.5)
    axs[0, 1].set_title('East position error, m')
    
    axs[0, 2].plot(time, down_pos_error, color=[0.45, 0, 0.9], linewidth=1.5)
    axs[0, 2].set_title('Down position error, m')
    
    axs[1, 0].plot(time, north_vel_error, color=[0.9, 0, 0.45], linewidth=1.5)
    axs[1, 0].set_title('North velocity error, m/s')
    
    axs[1, 1].plot(time, east_vel_error, color=[0.45, 0.9, 0], linewidth=1.5)
    axs[1, 1].set_title('East velocity error, m/s')
    
    axs[1, 2].plot(time, down_vel_error, color=[0, 0.45, 0.9], linewidth=1.5)
    axs[1, 2].set_title('Down velocity error, m/s')
    
    axs[2, 0].plot(time, roll_error, color=[0, 0.7, 0.7], linewidth=1.5)
    axs[2, 0].set_title('Attitude error about North, deg')
    axs[2, 0].set_xlabel('Time, s')
    
    axs[2, 1].plot(time, pitch_error, color=[0.7, 0, 0.7], linewidth=1.5)
    axs[2, 1].set_title('Attitude error about East, deg')
    axs[2, 1].set_xlabel('Time, s')
    
    axs[2, 2].plot(time, yaw_error, color=[0.7, 0.7, 0], linewidth=1.5)
    axs[2, 2].set_title('Heading error, deg')
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
