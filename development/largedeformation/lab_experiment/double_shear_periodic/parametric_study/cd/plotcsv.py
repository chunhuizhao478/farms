import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('code_cd_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data[:, 0], data[:, 1]/0.05, label='x', color='blue')
plt.show()