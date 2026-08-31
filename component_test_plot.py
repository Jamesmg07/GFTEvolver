import os
import re
import glob

import numpy as np
import matplotlib.pyplot as plt

Nx = 128
Ny = 128
Nz = 128

# Number of grid points in one density field
grid_size = Nx * Ny * Nz

data_dir = "25_08_0_0_tracking_test/Continual"
output_dir = data_dir

# ============================================================
# Which local density fields are enabled, and in what order.
# This must match the order local blocks are written in by
# Energy.cpp (energy, potential, gradient, kinetic, magnetic,
# electric) - only list the ones actually set to 1 in
# Energy.cfg, in their original relative order.
#
# Current Energy.cfg only has:
#   Output energy density? : 1
# with potential/gradient/kinetic/magnetic/electric density
# all set to 0, so only one block is written per file.
# ============================================================

enabled_densities = [
    "Energy density",
]

# ============================================================
# Find every local_quantities_<timestep>.dat file and sort by
# timestep (not alphabetically, since e.g. "1200" would
# otherwise sort before "300").
# ============================================================

filename_pattern = re.compile(r"local_quantities_(\d+)\.dat$")

file_paths = [
    path
    for path in glob.glob(os.path.join(data_dir, "local_quantities_*.dat"))
    if filename_pattern.search(os.path.basename(path))
]

file_paths.sort(
    key=lambda path: int(filename_pattern.search(os.path.basename(path)).group(1))
)

if not file_paths:
    raise FileNotFoundError(
        f"No local_quantities_*.dat files found in {data_dir}"
    )

# Centre of grid
x = Nx // 2
y = Ny // 2

z_values = np.arange(Nz)

# ============================================================
# Loop over every timestep's file and plot each enabled density
# ============================================================

for file_path in file_paths:

    timestep = int(
        filename_pattern.search(os.path.basename(file_path)).group(1)
    )

    data = np.loadtxt(file_path)

    densities = {}
    for block_iter, name in enumerate(enabled_densities):
        densities[name] = data[
            block_iter * grid_size : (block_iter + 1) * grid_size
        ].reshape((Nx, Ny, Nz))

    for name, density in densities.items():

        # --------------------------------------------------------
        # XZ slice at central y
        # --------------------------------------------------------

        slice_xz = density[:, y, :]

        plt.figure(figsize=(7, 6))

        plt.imshow(
            slice_xz,
            origin="lower",
            extent=[0, Nz - 1, 0, Nx - 1],
            aspect="equal"
        )

        plt.colorbar(label=name)

        plt.xlabel("z")
        plt.ylabel("x")
        plt.title(f"{name} at y = {y}, timestep {timestep}")

        plt.tight_layout()

        filename = name.lower().replace(" ", "_")

        plt.savefig(
            os.path.join(
                output_dir,
                f"{filename}_xz_timestep_{timestep}.png"
            ),
            dpi=300
        )

        plt.close()

        # --------------------------------------------------------
        # Density along central z line
        # --------------------------------------------------------

        density_centre_line = density[x, y, :]

        plt.figure(figsize=(7, 5))

        plt.plot(
            z_values,
            density_centre_line
        )

        plt.xlabel("z")
        plt.ylabel(name)

        plt.title(
            f"{name} along centre line "
            f"at (x, y) = ({x}, {y}), timestep {timestep}"
        )

        plt.tight_layout()

        plt.savefig(
            os.path.join(
                output_dir,
                f"{filename}_centre_line_timestep_{timestep}.png"
            ),
            dpi=300
        )

        plt.close()

print(f"Generated plots for {len(file_paths)} timesteps.")
