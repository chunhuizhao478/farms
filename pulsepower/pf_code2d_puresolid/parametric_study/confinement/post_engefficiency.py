import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools

total_energy = 20

DATASETS = {
    # Format: "Label": (path, linestyle, color)
    "Confinement: 1MPa": ("./case_1mpa/elasticity_csv.csv", "-", "b"),
    "Confinement: 5MPa": ("./case_5mpa/elasticity_csv.csv", "-", "r"),
    "Confinement: 10MPa": ("./case_10mpa/elasticity_csv.csv", "-", "g"),
}

# Column filtering options (choose one approach):
# Option 1: Hide specific columns (set column names to hide)
HIDDEN_COLUMNS = {"damping_work"}

# Option 2: Use only specific columns (if not empty, only these will be plotted)
# Leave empty to plot all columns (except hidden ones)
# Example: USE_COLUMNS = {"elastic_energy", "fracture_energy"}
USE_COLUMNS = {"dissipated_energy_total"}
# USE_COLUMNS = set()

# Color mode: "by_column" or "by_case"
# "by_column": Same column across cases gets same color (default)
# "by_case": Same case (dataset) gets same color across all columns
COLOR_MODE = "by_case"

time_column = "time"
t_max = 2e-5

# Pulse efficiency calculation settings
pulse_interval = 10e-6  # 10 microseconds in seconds
energy_per_pulse = 10   # J/pulse

# Create two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

default_colors = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
color_map_column = {}  # Map for "by_column" mode
color_map_case = {}  # Map for "by_case" mode

# Marker styles for different cases in scatter plot
markers = itertools.cycle(['o', 's', '^', 'D', 'v', '<', '>', 'p', '*', 'h'])
marker_map = {}

for label, (path, linestyle, user_color) in DATASETS.items():
    df = pd.read_csv(path, comment="/", skip_blank_lines=True)
    df = df[df[time_column] <= t_max].sort_values(time_column)
    time = df[time_column].to_numpy() * 1e6  # Convert to microseconds

    # Assign color and marker for this case
    if COLOR_MODE == "by_case":
        # Use user-specified color if in "by_case" mode
        color_map_case[label] = user_color

    if label not in marker_map:
        marker_map[label] = next(markers)

    for col in df.columns:
        if col == time_column or col in HIDDEN_COLUMNS:
            continue

        # If USE_COLUMNS is specified, only plot those columns
        if USE_COLUMNS and col not in USE_COLUMNS:
            continue

        # Determine color based on mode
        if COLOR_MODE == "by_case":
            plot_color = color_map_case[label]
        else:  # "by_column" mode (default)
            if col not in color_map_column:
                color_map_column[col] = next(default_colors)
            plot_color = color_map_column[col]

        y = df[col].to_numpy()
        time_s = df[time_column].to_numpy()  # Time in seconds for calculations
        mask = np.isfinite(time) & np.isfinite(y)

        # Plot 1: Time history
        ax1.plot(
            time[mask],
            y[mask],
            lw=1.0,
            linestyle=linestyle,
            color=plot_color,
            label=f"{label}",
        )

        # Calculate per-pulse efficiency for Plot 2
        pulse_times = []
        pulse_numbers = []
        pulse_efficiencies = []

        num_pulses = int(np.floor(t_max / pulse_interval))
        for i in range(1, num_pulses + 1):
            pulse_time = i * pulse_interval

            # Find dissipated energy at this pulse time
            idx = np.argmin(np.abs(time_s - pulse_time))
            energy_current = y[idx]

            # Calculate efficiency for this pulse
            if i == 1:
                energy_pulse = energy_current
            else:
                prev_pulse_time = (i - 1) * pulse_interval
                idx_prev = np.argmin(np.abs(time_s - prev_pulse_time))
                energy_prev = y[idx_prev]
                energy_pulse = energy_current - energy_prev

            efficiency = (energy_pulse / energy_per_pulse) * 100
            pulse_numbers.append(i)
            pulse_efficiencies.append(efficiency)

        # Plot 2: Scatter plot of per-pulse efficiency
        ax2.scatter(
            pulse_numbers,
            pulse_efficiencies,
            marker=marker_map[label],
            s=100,
            color=plot_color,
            label=f"{label}",
            edgecolors='black',
            linewidths=1.0,
        )

# Configure Plot 1 (Time history)
ax1.set_xlabel("time ($\mu$s)", fontsize=18)
ax1.set_ylabel("dissipated energy (J)", fontsize=18)
ax1.set_title("Pure Solid Dissipated Energy Time History", fontsize=20)
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend(loc="best", ncol=1, fontsize=12)
ax1.set_xlim(0, t_max * 1e6)
ax1.tick_params(axis="both", which="major", labelsize=14)

# Configure Plot 2 (Per-pulse efficiency)
ax2.set_xlabel("pulse number", fontsize=18)
ax2.set_ylabel("efficiency (%)", fontsize=18)
ax2.set_title("Per-Pulse Energy Efficiency", fontsize=20)
ax2.grid(True, ls=":", alpha=0.6)
ax2.legend(loc="best", ncol=1, fontsize=12)
ax2.tick_params(axis="both", which="major", labelsize=14)
ax2.set_xticks(range(1, num_pulses + 1))

plt.tight_layout()
plt.show()
