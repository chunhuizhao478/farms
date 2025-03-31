import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('code_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data[:, 0], data[:, 1]/0.05, label='x', color='blue')
plt.show()

#compute rate
rate = np.zeros((len(data)-1, 2))
for i in range(len(data)-1):
    rate[i, 0] = data[i+1, 0]
    rate[i, 1] = (data[i+1, 1] - data[i, 1]) / (data[i+1, 0] - data[i, 0])
plt.figure(figsize=(10, 6))
plt.plot(rate[:, 0], rate[:, 1], label='x', color='blue')
plt.xlabel('time')
plt.ylabel('rate')
plt.title('Rate of stress')