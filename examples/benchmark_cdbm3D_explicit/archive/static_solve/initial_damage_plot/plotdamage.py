import numpy as np
import matplotlib.pyplot as plt

# Parameters
e_damage = 0.7
e_sigma = 1000
x0, y0 = 0, 0  # Coordinates of the selected point

# Define the grid
x = np.linspace(-60000, 60000, 601)
y = np.linspace(-60000, 60000, 601)
X, Y = np.meshgrid(x, y)

# Calculate the distance from the selected point
r = np.sqrt((X - x0)**2 + (Y - y0)**2)

# Calculate the exponential decay
Z = e_damage * np.exp(-(r**2) / (e_sigma**2))

# Plotting
plt.figure()
plt.contourf(X, Y, Z, levels=50, cmap='viridis')
plt.colorbar(label='Damage Intensity')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('2D Exponential Decay from a Selected Point')
plt.show()

