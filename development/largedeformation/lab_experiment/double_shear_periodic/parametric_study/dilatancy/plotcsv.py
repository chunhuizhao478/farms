import numpy as np
import matplotlib.pyplot as plt

data0 = np.loadtxt('code_cg_csv_cg1em11.csv', delimiter=',', skiprows=1)
data1 = np.loadtxt('code_dilatancy_csv_1.csv', delimiter=',', skiprows=1)
data2 = np.loadtxt('code_dilatancy_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data0[:, 2], data0[:, 1]/0.05, label='x', color='blue')
plt.plot(data1[:, 2], data1[:, 1]/0.05, label='x', color='red')
plt.plot(data2[:, 2], data2[:, 1]/0.05, label='x', color='black')
plt.show()