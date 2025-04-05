import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('code_param_constraint_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(-data[:, 2]/0.095 * 100, data[:, 1]/(np.pi*0.027**2)/1e6, label='x', color='blue')
plt.xlabel('Axial strain (%)', fontsize=14)
plt.ylabel('Stress (MPa)', fontsize=14)
plt.title('Stress vs Axial Strain of Uniaxial Compression Test', fontsize=16)
plt.savefig('uniaxial_compression_test.png', dpi=300)
plt.show()