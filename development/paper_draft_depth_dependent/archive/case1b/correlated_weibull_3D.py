
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, weibull_min
from scipy.ndimage import gaussian_filter

# ---------------------------
# Step 0: Grid Setup
# ---------------------------
# Define a 3D spatial domain with physical coordinates
nx, ny, nz = 200, 20, 200
x = np.linspace(-20000, 20000, nx)
y = np.linspace(-1000, 1000, ny)
z = np.linspace(-20000, 0, nz)
xv, yv, zv = np.meshgrid(x, y, z, indexing='ij')  # not used in computation, just for reference

# ---------------------------
# Step 1: Generate Gaussian Random Field
# ---------------------------
# Start with standard Gaussian white noise (zero mean, unit variance)
# This field is uncorrelated (i.e., no spatial relationship between values)
np.random.seed(0)
gaussian_field = np.random.normal(0, 1, size=(nx, ny, nz))

# ---------------------------
# Step 2: Apply Spatial Correlation
# ---------------------------
# Use a Gaussian filter to smooth the field, introducing spatial correlation
# The sigma parameter controls the correlation length (in grid points)
correlation_length = 200
correlated_field = gaussian_filter(gaussian_field, sigma=correlation_length)

# ---------------------------
# Step 3: Gaussian to Uniform Mapping
# ---------------------------
# Use the CDF of the standard normal distribution to convert to a uniform field
# This maps values to [0, 1] while preserving spatial structure
uniform_field = norm.cdf(correlated_field)

# ---------------------------
# Step 4: Uniform to Weibull Mapping
# ---------------------------
# Apply the inverse CDF (percent point function) of the Weibull distribution
# This transforms the uniform values into Weibull-distributed values
# Parameters:
#   shape (k): controls the spread/steepness
#   scale (λ): controls the characteristic value
#   location (θ): minimum bound shift
weibull_shape = 2.0
weibull_scale = 0.1
weibull_location = 0.0
weibull_field = weibull_min.ppf(uniform_field, c=weibull_shape,
                                 scale=weibull_scale, loc=weibull_location)

# ---------------------------
# Step 5: Visualization
# ---------------------------
# Show a slice at the middle of the Y dimension using physical coordinates
y_index = ny // 2
plt.figure(figsize=(8, 6))
plt.pcolormesh(z, x, weibull_field[:, y_index, :], shading='auto', cmap='viridis')
plt.colorbar(label='Weibull Value')
plt.title(f'3D Correlated Weibull Field Slice (Y = {y[y_index]:.1f} m)')
plt.xlabel('Z [m]')
plt.ylabel('X [m]')
plt.tight_layout()
plt.show()