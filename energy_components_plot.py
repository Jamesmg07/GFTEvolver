import os

import numpy as np
import matplotlib.pyplot as plt

input_path = "Data/global_quantities.dat"
output_dir = os.path.dirname(input_path)

data = np.loadtxt(input_path)

# ============================================================
# global_quantities.dat column layout (one row per output):
#
# 0 = timestep
#
# followed by the enabled quantities, in this fixed order
# (only the enabled ones are written):
# 1 = Energy
# 2 = Potential energy
# 3 = Gradient energy
# 4 = Kinetic energy
# 5 = Magnetic energy
# 6 = Electric energy
#
# The config has all six enabled, output every 100 timesteps.
# ============================================================

timestep = data[:, 0]

energy = data[:, 1]
potential = data[:, 2]
gradient = data[:, 3]
kinetic = data[:, 4]
magnetic = data[:, 5]
electric = data[:, 6]

components = {
    "Energy": energy,
    "Potential energy": potential,
    "Gradient energy": gradient,
    "Kinetic energy": kinetic,
    "Magnetic energy": magnetic,
    "Electric energy": electric,
}

# ============================================================
# Combined plot: all energy components vs time
# ============================================================

plt.figure(figsize=(8, 6))

for name, values in components.items():
    plt.plot(timestep, values, label=name)

plt.xlabel("Timestep")
plt.ylabel("Energy")
plt.title("Energy components vs timestep")
plt.legend()
plt.tight_layout()

plt.savefig(
    os.path.join(output_dir, "energy_components_vs_time.png"),
    dpi=300
)

plt.close()

# ============================================================
# Individual plots, one per component (useful when components
# differ by orders of magnitude and don't compare well on a
# single shared axis)
# ============================================================

for name, values in components.items():

    plt.figure(figsize=(7, 5))

    plt.plot(timestep, values)

    plt.xlabel("Timestep")
    plt.ylabel(name)
    plt.title(f"{name} vs timestep")
    plt.tight_layout()

    filename = name.lower().replace(" ", "_")

    plt.savefig(
        os.path.join(output_dir, f"{filename}_vs_time.png"),
        dpi=300
    )

    plt.close()

print("All energy-component-vs-time plots generated.")
