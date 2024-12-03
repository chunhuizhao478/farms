import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.special import kv as modified_bessel_second_kind
from scipy.special import gamma as gamma_function
from scipy.interpolate import griddata
import pandas as pd
import netCDF4

# Grid parameters
nx, ny = 150, 80
x = np.linspace(-7500, 7500, nx)
y = np.linspace(-4000, 4000, ny)
x_grid, y_grid = np.meshgrid(x, y)
grid_points = np.column_stack((x_grid.ravel(), y_grid.ravel()))

# Von Kármán correlation parameters
correlation_length = 200
nu = 0.5

# Compute pairwise distances
distances = cdist(grid_points, grid_points, metric='euclidean')

# Von Kármán correlation function
def von_karman_correlation(r, l, nu):
    r_scaled = r / l
    with np.errstate(divide='ignore', invalid='ignore'):
        result = (2 ** (1 - nu)) / gamma_function(nu) * (r_scaled ** nu) * modified_bessel_second_kind(nu, r_scaled)
        result[np.isnan(result)] = 1.0
    return result

# Compute correlation matrix
correlation_matrix = von_karman_correlation(distances, correlation_length, nu)
np.fill_diagonal(correlation_matrix, 1.0)

# Generate random field
L = np.linalg.cholesky(correlation_matrix + 1e-6 * np.eye(len(grid_points)))
z = np.random.randn(len(grid_points))
von_karman_field = L @ z

# Reshape field to grid
von_karman_field_grid = von_karman_field.reshape(ny, nx)

# Shear modulus transformation
mean_shear_modulus = 10e9
std_dev_shear_modulus = 1e9
shear_modulus_field = mean_shear_modulus + std_dev_shear_modulus * von_karman_field_grid

# Parameters
file_path = './static_solve_out.e'

# Step 6: Map Values to Specific Points from NetCDF
# Get the data from the netCDF file
nc = netCDF4.Dataset(file_path)

# Flatten and prepare data for CSV
points = np.column_stack((x_grid.ravel(),y_grid.ravel()))
values = shear_modulus_field.ravel()

# Extract specific points from netCDF coordinates
x_coord = nc.variables['coordx'][:]
y_coord = nc.variables['coordy'][:]
target_points = np.column_stack((x_coord, y_coord))  # Target points for interpolation

# Define bounds
x_min, x_max = -7500, 7500
y_min, y_max = -4000, 4000

# Check if the points are within bounds
in_bounds = (x_coord >= x_min) & (x_coord <= x_max) & (y_coord >= y_min) & (y_coord <= y_max)

# Initialize mapped values with zeros
mapped_values = np.zeros_like(x_coord)

# Perform interpolation only for in-bound points
mapped_values[in_bounds] = griddata(points, values, target_points[in_bounds], method='linear')

# Combine the coordinates and mapped values into a single array
output_data = np.column_stack((x_coord, y_coord, mapped_values))

# Save the data to a CSV file
output_csv_file = "./mapped_vonkarman_field.csv"
np.savetxt(output_csv_file, output_data, delimiter=",", comments="")
print(f"Mapped data saved to {output_csv_file}")

# Plot the field
plt.figure(figsize=(10, 6))
plt.contourf(x, y, shear_modulus_field, levels=50, cmap='viridis')
plt.colorbar(label='Shear Modulus (Pa)')
plt.title('Von Kármán Random Distribution of Shear Modulus')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.show()