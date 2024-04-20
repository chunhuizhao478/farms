import numpy as np
import matplotlib.pyplot as plt

time1 = np.loadtxt("tdata_flat.txt")
# time2 = np.loadtxt("./outputs/list_timeseries0_2.txt")

pressure1 = np.loadtxt("pressure.txt")
# pressure2 = np.loadtxt("./outputs/list_pressure0_2.txt")

#plot
#global
plt.figure(figsize=(10,6))
plt.plot(time1[:]/(24*60*60),pressure1[:]/1e6,'b-')
# plt.plot(time2[:],pressure2[:]/1e6,'b-')
# plt.plot(timeseries[:85],ptr_y_valhist[:85],'k-',label = "velocity y")    
# plt.plot(timeseries[:],ptr_y_valhist,'r-',label = "velocity y")
plt.xlabel("time history (day)")
plt.ylabel("pressure (MPa)")
plt.title("Time History of Pressure at injection location with 10kg/s rate")
# plt.xlim([0,35])
# plt.ylim([0,18.0])
plt.savefig('img.png')

