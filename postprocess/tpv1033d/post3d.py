import numpy as np
import matplotlib.pyplot as plt

##loop time
num_begin = 0
num_end = 60
file_nums = np.linspace(num_begin,num_end,num_end+1)

time_begin = 0
time_end = 1.5
time = np.linspace(time_begin,time_end,num_end+1)

#data loc
loc = 502 - 1

#read data
##benchmark
data_uguca0 = np.loadtxt("./uguca_data/0strike7dot5dip.txt")
data_index_end = int(500)

data_uguca1 = np.loadtxt("./uguca_data/neg9strike7dot5dip.txt")
data_index_end = int(500) #DFM

data_uguca2 = np.loadtxt("./uguca_data/pos9strike7dot5dip.txt")
data_index_end = int(500) #DFM

timehist_strikeFaultMod = data_uguca0[:data_index_end+1,0]

sliprate_strikeFaultMod0 = data_uguca0[:data_index_end+1,2]
sliprate_strikeFaultMod1 = data_uguca1[:data_index_end+1,2]
sliprate_strikeFaultMod2 = data_uguca2[:data_index_end+1,2]

time_moose = np.loadtxt("./files/200m_929/time.txt")

sliprate_moose0 = np.loadtxt("./files/200m_929/sliprate_0strike0dot75dip.txt")
sliprate_moose1 = np.loadtxt("./files/200m_929/sliprate_neg9strike0dot75dip.txt")
sliprate_moose2 = np.loadtxt("./files/200m_929/sliprate_pos9strike0dot75dip.txt")

plt.figure()
plt.plot(time_moose,sliprate_moose0,'b--',label="moose-100m-3d")
plt.plot(timehist_strikeFaultMod,sliprate_strikeFaultMod0,'r-',label='uguca-100m-3d')
plt.title("sliprate along strike 0km dip 7.5km")
plt.xlabel("time (s)")
plt.ylabel("sliprate (m/s)")
plt.legend()

# plt.figure()
# plt.plot(time_moose,sliprate_moose1,'b--',label="moose")
# plt.plot(timehist_strikeFaultMod,sliprate_strikeFaultMod1,'r-',label='uguca')
# plt.title("sliprate along strike -9km dip 7.5km")
# plt.xlabel("time (s)")
# plt.ylabel("sliprate (m/s)")
# plt.legend()

# plt.figure()
# plt.plot(time_moose,sliprate_moose2,'b--',label="moose")
# plt.plot(timehist_strikeFaultMod,sliprate_strikeFaultMod2,'r-',label='uguca')
# plt.title("sliprate along strike 9km dip 7.5km")
# plt.xlabel("time (s)")
# plt.ylabel("sliprate (m/s)")
# plt.legend()

#
plt.show()