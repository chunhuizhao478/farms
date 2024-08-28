import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Domain dimensions
x_size = 120
y_size = 60
z_size = 120

# Create coordinate grids
x = np.linspace(-x_size/2, x_size/2, x_size+1)
y = np.linspace(0, -y_size, y_size+1)
z = np.linspace(-z_size/2, z_size/2, z_size+1)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Gaussian parameters
peak_value = 0.7
y_center = -7.5
sigma_x = x_size / 10  # Controls spread in x-direction
sigma_y = y_size / 20  # Controls spread in y-direction
sigma_z = z_size / 10  # Controls spread in z-direction

# 3D Gaussian function
gaussian_3d = peak_value * np.exp(-(X**2 / (2 * sigma_x**2) +
                                    (Y-y_center)**2 / (2 * sigma_y**2) +
                                    Z**2 / (2 * sigma_z**2)))

# Fixed X and Y values (center in X, center in Y where Gaussian peaks)
x_fixed = 0
y_fixed = -7.5

# Find the index closest to x_fixed and y_fixed
x_index = np.argmin(np.abs(x - x_fixed))
y_index = np.argmin(np.abs(y - y_fixed))

# Extract the Gaussian values along the Z direction at the fixed X and Y
gaussian_line = gaussian_3d[x_index, y_index, :]

# Plot the Gaussian values along Z
plt.figure(figsize=(10, 6))
plt.plot(z, gaussian_line, label=f'Gaussian at X={x_fixed}, Y={y_fixed}', color='blue')
plt.xlabel('Z')
plt.ylabel('Gaussian Value')
plt.title(f'Gaussian Distribution along Z at X={x_fixed}, Y={y_fixed}')
plt.grid(True)
plt.legend()
plt.show()
