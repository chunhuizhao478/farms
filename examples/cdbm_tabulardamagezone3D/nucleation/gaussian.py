import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt

# Define the grid size and parameters
L = 4  # Length of the square box
grid_size = 10  # Number of points in each dimension
A = 1.0  # Peak value at the center
sigma = 0.8  # Spread of the Gaussian

# Create the grid
x = np.linspace(-L/2, L/2, grid_size)
y = np.linspace(-L/2, L/2, grid_size)
X, Y = np.meshgrid(x, y)

# Calculate the distance from the center
D = np.sqrt(X**2 + Y**2)

# Apply the Gaussian function
V = A * np.exp(-D**2 / (2 * sigma**2))

# Theoretical min and max values of V
V_max = A
d_max = (L * np.sqrt(2)) / 2
V_min = A * np.exp(-(d_max**2) / (2 * sigma**2))

# Normalize V to range [0, 1] using theoretical min and max
V_normalized = (V - V_min) / (V_max - V_min)

# Scale V to range [0.7, 0.8]
V_scaled = V_normalized * (0.8 - 0.7) + 0.7

# Plot the result
plt.imshow(V_scaled, extent=[-L/2, L/2, -L/2, L/2], origin='lower')
plt.colorbar(label='Value')
plt.title('Values Highest at Center and Gradually Decreasing to Boundary')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

