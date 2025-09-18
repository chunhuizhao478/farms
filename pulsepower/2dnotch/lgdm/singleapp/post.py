import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("elasticity_csv.csv", delimiter=",", skiprows=1)

rate = 4e-3

plt.figure(figsize=(10, 6))
plt.plot(data[:, 0] * rate, data[:, 1], 'r-')
plt.xlabel("Displacement (m)")
plt.ylabel("Force (N)")
plt.title("Reaction Force vs Displacement")
plt.grid()
plt.show()