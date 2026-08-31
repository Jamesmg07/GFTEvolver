import os

import numpy as np
import matplotlib.pyplot as plt

input_path = "25_08_0_0_tracking_test/monopole_tracking.dat"
output_dir = os.path.dirname(input_path)

data = np.loadtxt(input_path)

# ============================================================
# monopole_tracking.dat column layout (one row per output):
#
# 0 = timestep
# 1,2,3,4 = x, y, z, value   for monopole 1 (largest local max)
# 5,6,7,8 = x, y, z, value   for monopole 2 (second local max,
#                              at least "Minimum physical
#                              monopole separation" away from
#                              monopole 1)
#
# The C++ analyser writes first.x*dx, first.y*dy, first.z*dz
# (see GaugeCondition/Energy analyser source), so the x, y, z
# columns are already physical distances, not grid indices -
# no further scaling is needed here. "value" is the local
# energy-density value at that point (an intensive quantity),
# so it's physical as-is too.
#
# Rows where a monopole wasn't found have "nan" in the
# corresponding columns.
# ============================================================

timestep = data[:, 0]

x1, y1, z1, value1 = data[:, 1], data[:, 2], data[:, 3], data[:, 4]
x2, y2, z2, value2 = data[:, 5], data[:, 6], data[:, 7], data[:, 8]

# Minimum physical monopole separation enforced by the analyser
# (from the config file). Used as a reference line only.
minimum_separation = 3.0

# --------------------------------------------------------
# Physical separation between the two tracked monopoles.
# NaNs propagate naturally, so time steps where either
# monopole wasn't found simply show up as gaps in the line.
# --------------------------------------------------------

separation = np.sqrt(
    (x2 - x1) ** 2
    + (y2 - y1) ** 2
    + (z2 - z1) ** 2
)

# ============================================================
# Plot separation vs time
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(
    timestep,
    separation,
    marker="o",
    markersize=2,
    linewidth=1,
    label="Monopole separation"
)

plt.axhline(
    minimum_separation,
    color="grey",
    linestyle="--",
    linewidth=1,
    label=f"Minimum enforced separation ({minimum_separation})"
)

plt.xlabel("Timestep")
plt.ylabel("Physical separation")
plt.title("Monopole separation vs timestep")
plt.legend()
plt.tight_layout()

plt.savefig(
    os.path.join(output_dir, "monopole_separation_vs_time.png"),
    dpi=300
)

plt.close()

# ============================================================
# Bonus: monopole peak values vs time, since they're already
# loaded and often useful alongside separation.
# ============================================================

plt.figure(figsize=(7, 5))

plt.plot(timestep, value1, label="Monopole 1 peak value")
plt.plot(timestep, value2, label="Monopole 2 peak value")

plt.xlabel("Timestep")
plt.ylabel("Peak local energy density")
plt.title("Monopole peak values vs timestep")
plt.legend()
plt.tight_layout()

plt.savefig(
    os.path.join(output_dir, "monopole_peak_values_vs_time.png"),
    dpi=300
)

plt.close()

print("Monopole separation plots generated.")
