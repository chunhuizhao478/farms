import numpy as np
import matplotlib.pyplot as plt

# Physical constants
density_fluid = 1000      # kg/m^3
rho = 2670                # kg/m^3
g = 9.8                   # m/s^2

# Elastic constants
mu = 32.04e9              # shear modulus (Pa)
lmbda = 32.04e9           # Lame's first parameter (Pa)

# Depth array (m)
depths = np.linspace(0, 15000, 300)  # from surface to 15 km

# Pore pressure and vertical stress
Pf = density_fluid * g * depths
sigma_zz = -rho * g * depths

# Coefficients for horizontal and shear stresses
b_xx = 3.5
b_yy = 1.0
b_xy = -0.6

# Piecewise definitions
mask = depths <= 15600
sigma_xx = np.where(mask, b_xx * (sigma_zz + Pf) - Pf, sigma_zz)
sigma_yy = np.where(mask, b_yy * (sigma_zz + Pf) - Pf, sigma_zz)
sigma_xy = np.where(mask, b_xy * (sigma_zz + Pf), 0.0)

# cohesion
mask = depths <= 4000
c = np.where(mask, 0.3e6 + (0.000675e6) * (4000 - depths), 0.3e6)  # cohesion in Pa

# shear strength
mu_s = 0.677
static_shear_strength = c + abs( mu_s * (sigma_yy + Pf) )
mu_d = 0.525
residual_shear_strength = c + abs( mu_d * (sigma_yy + Pf) )

# ------------------------
# Plot Stress Components
# ------------------------
plt.figure(figsize=(6, 8))
plt.plot(Pf/1e6, depths/1e3, label=r'$P_f$')
plt.plot(abs(sigma_zz)/1e6, depths/1e3, label=r'$\sigma_{zz}$') #when plotting, we use positive values for compression
plt.plot(abs(sigma_xx)/1e6, depths/1e3, label=r'$\sigma_{xx}$')
plt.plot(abs(sigma_yy)/1e6, depths/1e3, label=r'$\sigma_{yy}$')
plt.plot(sigma_xy/1e6, depths/1e3, label=r'$\sigma_{xy}$')
plt.plot(static_shear_strength/1e6, depths/1e3, label='Static Shear Strength', linestyle='--', color='orange')
plt.plot(residual_shear_strength/1e6, depths/1e3, label='Residual Shear Strength', linestyle='--', color='red')
plt.gca().invert_yaxis()
plt.ylabel('Depth (km)')
plt.xlabel('Stress (MPa)')
plt.title('Stress Components vs Depth')
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.savefig('stress_components_vs_depth.png', dpi=300)
plt.show()

# ------------------------
# Compute Strains
# ------------------------
# Assemble stress tensor components
stress = np.zeros((len(depths), 3, 3))
stress[:, 0, 0] = sigma_xx
stress[:, 1, 1] = sigma_yy
stress[:, 2, 2] = sigma_zz
stress[:, 0, 1] = sigma_xy
stress[:, 1, 0] = sigma_xy

# Inverse Hooke's law: epsilon_ij = 1/(2μ) σ_ij - λ/(2μ(3λ+2μ)) σ_kk δ_ij
trace_coeff = -lmbda / (2 * mu * (3 * lmbda + 2 * mu))
strain = np.zeros_like(stress)
enum = np.trace(stress, axis1=1, axis2=2)  # vector of σ_kk
for i in range(3):
    for j in range(3):
        strain[:, i, j] = stress[:, i, j] / (2 * mu)
for k in range(3):
    strain[:, k, k] += trace_coeff * enum

# Compute invariants: I1, I2, xi
I1 = np.trace(strain, axis1=1, axis2=2)
I2 = np.einsum('nij,nij->n', strain, strain)
xi = I1 / np.sqrt(I2)

# ------------------------
# Plot Strain Invariants
# ------------------------
plt.figure(figsize=(6, 8))
plt.plot(I1, depths, label=r'$I_1$')
plt.plot(I2, depths, label=r'$I_2$')
plt.plot(xi, depths, label=r'$\xi$')
plt.gca().invert_yaxis()
plt.ylabel('Depth (m)')
plt.xlabel('Invariant values')
plt.title('Strain Invariants vs Depth')
plt.legend(loc='best')
plt.grid(True)
plt.tight_layout()
plt.savefig('strain_invariants_vs_depth.png', dpi=300)
plt.show()