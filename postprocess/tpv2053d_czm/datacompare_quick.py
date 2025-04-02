import numpy as np
import matplotlib.pyplot as plt

benchmark_code = "TPV205"

benchmark_label = "benchmark-eqsim-200m"
farms_label = "farms-200m"

xcoord_i = 0
zcoord_i = -7500

xcoord_i = xcoord_i / 1e3
zcoord_i = zcoord_i / 1e3

benchmark_path = "./benchmark_data/eqsim_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
farms_path = "./farms_data_manual/czm_strike0.0_dip-7.5.txt"

benchmark = np.loadtxt(benchmark_path)
farms = np.loadtxt(farms_path, skiprows=1, delimiter=",")

## slip
plt.figure()
plt.plot(farms[:,0],-farms[:,1],'g-',label=farms_label)
plt.plot(benchmark[:,0],benchmark[:,1],'r-',label=benchmark_label)
plt.title(benchmark_code+" slip time history at strike "+str(xcoord_i)+"km and at dip "+str(zcoord_i)+"km ")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("slip (m)")
plt.show()

## slip rate
plt.figure()
plt.plot(farms[:,0],-farms[:,2],'g-',label=farms_label)
plt.plot(benchmark[:,0],benchmark[:,2],'r-',label=benchmark_label)
plt.title(benchmark_code+" slip rate time history at strike "+str(xcoord_i)+"km and at dip "+str(zcoord_i)+"km ")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("slip rate (m/s)")
plt.show()   