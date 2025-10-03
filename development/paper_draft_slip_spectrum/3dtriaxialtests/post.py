import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('./parametric_results/sigma_000MPa_main.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data[:, 0], data[:, 1]/(np.pi*0.027*0.027), "r-*")
plt.show()