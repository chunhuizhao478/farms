import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("elasticity_out.csv", skiprows=1, delimiter=",")

plt.figure(figsize=(10, 6))
plt.plot(data[:, 0] * 1e-7 * 1e3, data[:, 1] / 1e6, label="QUAD Elem", color='blue')
plt.xlabel("Displacement (mm)")
plt.ylabel("Stress (MPa)")
plt.title("Uniaxial Tension Test")
plt.show()