import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('dynamic_solve_main_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data[:, 0], data[:, 1]/0.01, "r-*")
plt.show()