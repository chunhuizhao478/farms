import numpy as np
import matplotlib.pyplot as plt

data_1 = np.loadtxt('code_velstep_csv_cg_1em11_m2_1.csv', delimiter=',', skiprows=1)
data_2 = np.loadtxt('code_velstep_csv_cg_1em11_m2_1.05.csv', delimiter=',', skiprows=1)
data_3 = np.loadtxt('code_velstep_csv_cg_1em11_m2_1.1.csv', delimiter=',', skiprows=1)
data_4 = np.loadtxt('code_velstep_csv_cg_1em11_m2_1.2.csv', delimiter=',', skiprows=1)
data_5 = np.loadtxt('code_velstep_csv_cg_1em12_m2_1.2.csv', delimiter=',', skiprows=1)
data_6 = np.loadtxt('code_velstep_csv.csv', delimiter=',', skiprows=1)
plt.figure(figsize=(10, 6))
# plt.plot(data_1[:, 2], data_1[:, 1]/0.05, label='x', color='blue')
plt.plot(data_2[:, 2], data_2[:, 1]/0.05, label='x', color='black')
plt.plot(data_3[:, 2], data_3[:, 1]/0.05, label='x', color='red')
plt.plot(data_4[:, 2], data_4[:, 1]/0.05, label='x', color='green')
plt.plot(data_5[:, 2], data_5[:, 1]/0.05, label='x', color='green')
plt.plot(data_6[:, 2], data_6[:, 1]/0.05, label='x', color='purple')
# plt.axhline(3.04e7, color='black', lw=1)
plt.xlim([0, 0.0004])
plt.ylim([0, 4e7])
plt.show()

# #compute rate
# rate = np.zeros((len(data)-1, 2))
# for i in range(len(data)-1):
#     rate[i, 0] = data[i+1, 2]
#     rate[i, 1] = (data[i+1, 1] - data[i, 1]) / 0.05
# plt.figure(figsize=(10, 6))
# plt.plot(rate[:, 0], rate[:, 1]/1e6, label='x', color='blue')
# # plt.axhline(0, color='black', lw=1)
# plt.xlabel('time')
# plt.ylabel('rate')
# plt.title('change of stress')
# plt.ylim([-0.01, 0.01])
# plt.show()