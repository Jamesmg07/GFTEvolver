import numpy as np

Nx = 256
Ny = 256
Nz = 256

data = np.loadtxt(
    "Data/Continual/local_quantities_0.dat"
)

energy = data.reshape((Nz, Ny, Nx))

y = Ny // 2
x = Nx // 2

slice_xz = energy[:, y, :]

import matplotlib.pyplot as plt

plt.figure()

plt.imshow(
    slice_xz,
    origin="lower",
    extent=[0, Nx, 0, Nz],
    aspect="equal"
)

plt.colorbar(label="Energy density")

plt.xlabel("z")
plt.ylabel("x")
plt.title(f"Energy density at y = {y}")

plt.show()

#plot along z centre line, energy density vs z

z_values = np.arange(Nz)


y = Ny // 2
x = Nx // 2

energy_centre_line = energy[x, y, :]

plt.figure()
plt.plot(z_values, energy_centre_line)
plt.xlabel("z")
plt.ylabel("Energy density")
plt.title(f"Energy density along centre line at (x, y) = ({x}, {y})")
plt.show()