import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('elasticity_csv.csv', delimiter=',', skiprows=1)

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(data[:, 0], data[:, 1], label='dissipated_energy_total')
plt.plot(data[:, 0], data[:, 2], label='full_energy')
plt.plot(data[:, 0], data[:, 3], label='full_input_energy')
plt.plot(data[:, 0], data[:, 4], label='solid_elastic_energy_total')
plt.plot(data[:, 0], data[:, 5], label='solid_kinetic_energy_total')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.title('Elasticity Data')
plt.legend()
plt.grid()
plt.show()