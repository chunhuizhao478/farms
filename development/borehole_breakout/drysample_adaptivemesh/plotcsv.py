import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.loadtxt('code_drysample_smallstrain_cd50_csv.csv', delimiter=',', skiprows=1)

# Compute strain (dimensionless) and stress (in MPa)
# Note: The strain is computed as -displacement/initial height.
strain = -data[:, 2] / 0.105  # dimensionless strain
stress = -data[:, 1] / (np.pi * 0.027**2) / 1e6  # stress in MPa

# Plot settings: converting strain to percentage for plotting
strain_percent = strain * 100

# Define a mask for the elastic (linear) region. Here we assume that strains up to 1% (0.01 in dimensionless form) are elastic.
mask = strain <= 0.01

# Perform a linear regression on the elastic region
coeff = np.polyfit(strain[mask], stress[mask], 1)
slope = coeff[0]      # This is the slope (Young's modulus in MPa)
intercept = coeff[1]  # Intercept of the fit

print("Computed slope (Young's modulus):", slope, "MPa")

# Calculate the fitted stress values over the full strain range for plotting
fitted_stress = np.polyval(coeff, strain)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(strain_percent, stress, label='Experimental Data', color='blue')
# plt.plot(strain_percent, fitted_stress, '--', label=f'Linear Fit (slope = {slope:.2f} MPa)', color='red')
plt.xlabel('Axial Strain (%)', fontsize=14)
plt.ylabel('Stress (MPa)', fontsize=14)
plt.title('Stress vs Axial Strain of Uniaxial Compression Test', fontsize=16)
plt.legend()
plt.savefig('uniaxial_compression_test.png', dpi=300)
plt.show()