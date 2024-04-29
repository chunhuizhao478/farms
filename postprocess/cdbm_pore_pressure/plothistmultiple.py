import numpy as np
import matplotlib.pyplot as plt

time1 = np.loadtxt("tdata_flat.txt")
# time2 = np.loadtxt("./outputs/list_timeseries0_2.txt")

pressure1 = np.loadtxt("pressure3em12.txt")
pressure2 = np.loadtxt("pressure3em13.txt")
pressure3 = np.loadtxt("pressure3em14.txt")
pressure4 = np.loadtxt("pressure3em15.txt")
pressure5 = np.loadtxt("pressure3em16.txt")

#plot
#global
plt.figure(figsize=(10,6))
plt.loglog(time1[:]/(24*60*60),pressure1[:]/1e6,'b-',label = "k = 3em12")
plt.loglog(time1[:]/(24*60*60),pressure2[:]/1e6,'k-',label = "k = 3em13")
plt.loglog(time1[:]/(24*60*60),pressure3[:]/1e6,'r-',label = "k = 3em14")  
plt.loglog(time1[:]/(24*60*60),pressure4[:]/1e6,'m-',label = "k = 3em15")
plt.loglog(time1[:]/(24*60*60),pressure5[:]/1e6,'g-',label = "k = 3em16")
plt.xlabel("time history (day)")
plt.ylabel("pressure (MPa)")
plt.title("Time History of Pressure at injection location with 10kg/s rate")
plt.legend()
# plt.xlim([0,35])
# plt.ylim([0,18.0])
plt.savefig('img.png')

