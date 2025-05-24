import numpy as np
import matplotlib.pyplot as plt

data_1 = np.loadtxt('code_velstep_v3em6_csv.csv', delimiter=',', skiprows=1)
data_2 = np.loadtxt('code_velstep_v5em6_csv.csv', delimiter=',', skiprows=1)
data_3 = np.loadtxt('code_velstep_v7em6_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data_1[:, 2], data_1[:, 1]/0.05, label='vel = 3e-6 m/s', color='blue')
plt.plot(data_2[:, 2], data_2[:, 1]/0.05, label='vel = 5e-6 m/s', color='black')
plt.plot(data_3[:, 2], data_3[:, 1]/0.05, label='vel = 7e-6 m/s', color='red')
# plt.axhline(3.04e7, color='black', lw=1)
plt.legend()
plt.title("Shear Stress-strain curve", fontsize=18)
plt.xlabel("Shear strain", fontsize=16)
plt.ylabel("Shear stress (MPa)", fontsize=16)
plt.xlim(0, 0.0008)
plt.show()
plt.savefig("velstep.png")