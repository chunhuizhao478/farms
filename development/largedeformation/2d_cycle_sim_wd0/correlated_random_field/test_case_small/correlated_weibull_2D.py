import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import scipy.special as sp_gamma
from scipy.spatial.distance import cdist
from scipy.interpolate import griddata
import netCDF4

# Parameters
file_path = './static_solve_out.e'

# Step 1: Define the Spatial Grid
nx = 150  # Number of grid points in each dimension
ny = 150
x = np.linspace(-7500, 7500, nx)
y = np.linspace(-7500, 7500, ny)
x_grid, y_grid = np.meshgrid(x, y)
grid_points = np.column_stack((x_grid.ravel(), y_grid.ravel()))  # Flatten grid into points

# Step 2: Compute the Covariance Matrix
correlation_length = 200  # Define correlation length
num_points = grid_points.shape[0]

# Compute pairwise distances efficiently
distances = cdist(grid_points, grid_points, metric='euclidean')  # Pairwise distances

# Apply the exponential decay function
covariance_matrix = np.exp(-distances / correlation_length)

# Step 3: Generate a Gaussian Random Field
# Add a small value to the diagonal for numerical stability
covariance_matrix += 1e-6 * np.eye(num_points)

# Cholesky decomposition
L = np.linalg.cholesky(covariance_matrix)

# Generate uncorrelated standard normal random variables
z = np.random.randn(num_points)

# Generate the correlated Gaussian random field
gaussian_field = L @ z
gaussian_field_grid = gaussian_field.reshape(ny, nx)  # Reshape to 2D grid

# Step 4: Transform to Weibull Distribution
# Mean friction angle (in degrees)
mu = 46.8

# Typical shape parameters for Weibull distribution (small k -> larger spread)
shape_param = 12

# Compute the scale parameter lambda for each shape parameter k
scale_param = mu / sp_gamma.gamma(1 + 1/shape_param)

print(f'Shape parameter k: {shape_param}')
print(f'Scale parameter λ: {scale_param}')

# Convert Gaussian values to uniform values using the standard normal CDF
uniform_values = stats.norm.cdf(gaussian_field)

# Apply the inverse Weibull CDF (quantile function)
weibull_field = scale_param * (-np.log(1 - uniform_values))**(1 / shape_param)
weibull_field_grid = weibull_field.reshape(ny, nx)  # Reshape to 2D grid

# Transform to xi_o
xi_o_field_grid = -np.sqrt(2) / np.sqrt(1 + (1 + 1) * (1 + 1) *
                                        np.sin(weibull_field_grid * np.pi / 180) *
                                        np.sin(weibull_field_grid * np.pi / 180))

# Step 5: Visualize the Results
# Plot Gaussian random field
plt.figure()
plt.contourf(x_grid, y_grid, gaussian_field_grid, cmap='viridis')
plt.title('Correlated Gaussian Random Field')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Amplitude')
plt.show()

# Plot Weibull random field
plt.figure()
plt.contourf(x_grid, y_grid, weibull_field_grid, cmap='viridis')
plt.title(r'Correlated Weibull Random Field $\text{Friction Angle \Phi}$')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Value')
plt.show()

# Plot xi_o field
plt.figure()
plt.contourf(x_grid, y_grid, xi_o_field_grid, cmap='viridis')
plt.title(r'Correlated Weibull Random Field $\text{\xi_o}$')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(label='Value')
plt.show()

# Step 6: Map Values to Specific Points from NetCDF
# Get the data from the netCDF file
nc = netCDF4.Dataset(file_path)

# Flatten the grid and Weibull field for interpolation
points = np.column_stack((x_grid.ravel(), y_grid.ravel()))  # All grid points
values = xi_o_field_grid.ravel()  # Corresponding values

# Extract specific points from netCDF coordinates
x_coord = nc.variables['coordx'][:]
y_coord = nc.variables['coordy'][:]
target_points = np.column_stack((x_coord, y_coord))  # Target points for interpolation

# Define bounds
x_min, x_max = -7500, 7500
y_min, y_max = -7500, 7500

# Check if the points are within bounds
in_bounds = (x_coord >= x_min) & (x_coord <= x_max) & (y_coord >= y_min) & (y_coord <= y_max)

# Initialize mapped values with zeros
mapped_values = np.zeros_like(x_coord)

# Perform interpolation only for in-bound points
mapped_values[in_bounds] = griddata(points, values, target_points[in_bounds], method='linear')

# Combine the coordinates and mapped values into a single array
output_data = np.column_stack((x_coord, y_coord, mapped_values))

# Save the data to a CSV file
output_csv_file = "./mapped_weibull_field.csv"
np.savetxt(output_csv_file, output_data, delimiter=",", comments="")
print(f"Mapped data saved to {output_csv_file}")