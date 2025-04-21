import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('code_velstep_wstatevariable_vw_csv.csv', delimiter=',', skiprows=1)
data2 = np.loadtxt('code_velstep_wstatevariable_vw_2_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(data[:, 2], data[:, 1]/0.05, 'b-')
plt.plot(data2[:, 2], data2[:, 1]/0.05, 'r-')
# plt.axhline(3.04e7, color='black', lw=1)
plt.show()

#compute rate
# rate = np.zeros((len(data)-1, 2))
# for i in range(len(data)-1):
#     rate[i, 0] = data[i+1, 2]
#     rate[i, 1] = (data[i+1, 1] - data[i, 1]) / 0.05
# plt.figure(figsize=(10, 6))
# # plt.plot(rate[:, 0], rate[:, 1]/1e6, label='x', color='blue')
# plt.axhline(0, color='black', lw=1)
# plt.xlabel('time')
# plt.ylabel('rate')
# plt.title('change of stress')
# plt.ylim([-0.01, 0.01])
# plt.show()