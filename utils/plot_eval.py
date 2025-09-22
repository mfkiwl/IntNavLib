import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

def plot_eval(csv_file):
    # Load data
    data = pd.read_csv(csv_file)
    time = data.iloc[:, 0].to_numpy()

    # ------------------------------------------------------------------
    # Utility function: plot curve with ± sigma bounds
    # ------------------------------------------------------------------
    def plot_with_bounds(ax, time, value, sigma, color, title, xlabel=None, centered=True):
        ax.plot(time, value, color=color, linewidth=1.5)
        if centered:
            # error plot (centered at 0)
            ax.plot(time, sigma, color=color, linestyle="dashed", linewidth=1.0)
            ax.plot(time, -sigma, color=color, linestyle="dashed", linewidth=1.0)
        else:
            # bias/estimate plot (centered at estimate)
            ax.plot(time, value + sigma, color=color, linestyle="dashed", linewidth=1.0)
            ax.plot(time, value - sigma, color=color, linestyle="dashed", linewidth=1.0)
        ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)

    # ------------------------------------------------------------------
    # Fixed block 1: navigation errors (position, velocity, attitude)
    # ------------------------------------------------------------------
    errors = [
        ("X position error, m",  1, 10, [0.9, 0.45, 0]),
        ("Y position error, m",  2, 11, [0, 0.9, 0.45]),
        ("Z position error, m",  3, 12, [0.45, 0, 0.9]),
        ("X velocity error, m/s",4, 13, [0.9, 0, 0.45]),
        ("Y velocity error, m/s",5, 14, [0.45, 0.9, 0]),
        ("Z velocity error, m/s",6, 15, [0, 0.45, 0.9]),
        ("X rotation error, deg",7, 16, [0, 0.7, 0.7]),
        ("Y rotation error, deg",8, 17, [0.7, 0, 0.7]),
        ("Z rotation error, deg",9, 18, [0.7, 0.7, 0]),
    ]

    fig, axs = plt.subplots(3, 3, figsize=(15, 10))
    fig.suptitle("Navigation Errors", fontsize=16)
    fig.canvas.manager.set_window_title("Navigation Errors")
    fig.tight_layout(pad=5.0)

    for i, (title, err_idx, sigma_idx, color) in enumerate(errors):
        row, col = divmod(i, 3)
        error = data.iloc[:, err_idx].to_numpy()
        sigma = data.iloc[:, sigma_idx].to_numpy() * 3
        xlabel = "Time, s" if row == 2 else None
        plot_with_bounds(axs[row, col], time, error, sigma, color, title, xlabel, centered=True)

    # ------------------------------------------------------------------
    # Fixed block 2: bias/clock estimates + sigmas
    # ------------------------------------------------------------------
    # order: 3 acc bias, 3 gyro bias, clock offset, clock drift
    # indices: 19–26 (estimates), 27–34 (sigmas)
    biases = [
        ("Accel bias X, m/s^2", 19, 27, [0.8, 0.2, 0.2]),
        ("Accel bias Y, m/s^2", 20, 28, [0.2, 0.8, 0.2]),
        ("Accel bias Z, m/s^2", 21, 29, [0.2, 0.2, 0.8]),
        ("Gyro bias X, rad/s",  22, 30, [0.8, 0.5, 0.2]),
        ("Gyro bias Y, rad/s",  23, 31, [0.5, 0.8, 0.2]),
        ("Gyro bias Z, rad/s",  24, 32, [0.2, 0.5, 0.8]),
        ("Clock offset, m", 25, 33, [0.6, 0.2, 0.6]),
        ("Clock drift, m/s",  26, 34, [0.2, 0.6, 0.6]),
    ]

    fig_bias, axs_bias = plt.subplots(4, 2, figsize=(12, 12))
    axs_bias = axs_bias.flatten()
    fig_bias.canvas.manager.set_window_title("IMU and Clock Bias Estimates")
    fig_bias.suptitle("IMU and Clock Bias Estimates", fontsize=16)
    fig_bias.tight_layout(pad=5.0)

    for i, (title, est_idx, sigma_idx, color) in enumerate(biases):
        estimate = data.iloc[:, est_idx].to_numpy()
        sigma = data.iloc[:, sigma_idx].to_numpy() * 3
        xlabel = "Time, s" if i >= 6 else None
        plot_with_bounds(axs_bias[i], time, estimate, sigma, color, title, xlabel, centered=False)

    # ------------------------------------------------------------------
    # Block 3: dynamic innovations
    # ------------------------------------------------------------------
    extra_cols = data.columns[35:]
    num_extra = len(extra_cols)

    if num_extra > 0:
        assert num_extra % 2 == 0, "Extra columns must come in innovation/sigma pairs"
        n_pairs = num_extra // 2
        ncols = 2
        nrows = math.ceil(n_pairs / ncols)

        fig_innov, axs_innov = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows))
        axs_innov = np.array(axs_innov).reshape(-1)
        fig_innov.canvas.manager.set_window_title("Kalman Filter Innovations")
        fig_innov.suptitle("Kalman Filter Innovations", fontsize=16)
        fig_innov.tight_layout(pad=5.0)

        for i in range(n_pairs):
            innov = data.iloc[:, 35 + 2*i].to_numpy()
            sigma = data.iloc[:, 35 + 2*i + 1].to_numpy() * 3
            ax = axs_innov[i]
            title = f"Innovation {i+1}"
            plot_with_bounds(ax, time, innov, sigma, [0.2, 0.2, 0.8], title, xlabel="Time, s", centered=True)

        # Hide unused axes if any
        for j in range(n_pairs, len(axs_innov)):
            fig_innov.delaxes(axs_innov[j])

    # Show all figures
    plt.show()


def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plot_eval.py <csv_file>")
        return
    plot_eval(sys.argv[1])

if __name__ == "__main__":
    main()
