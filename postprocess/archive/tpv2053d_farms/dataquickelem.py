import numpy as np
import matplotlib.pyplot as plt

#read benchmark data
benchmark_label = "benchmark-eqsim-200m"
benchmark_code = "tpv205"

#read farms data
farms_label = "farms-200m"

xcoord_i = 0.0;
zcoord_i = 0.0;

#time farms 0.1s interval
time = np.linspace(0,12.0,121)
datalen = 97

#file path
benchmark_path = "./benchmark_data/eqsim_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
farms_slip_data_path = "/Users/andyz/projects/farms_benchmark/postprocess/tpv2053d/zdir/jumpz_strike0dip0.txt"
farms_traction_data_path = "/Users/andyz/projects/farms_benchmark/postprocess/tpv2053d/zdir/tractionz_strike0dip0.txt"

#
benchmark = np.loadtxt(benchmark_path)
farms_slip_z = np.loadtxt(farms_slip_data_path,skiprows=1)
farms_traction_z = np.loadtxt(farms_traction_data_path,skiprows=1)

## slip z
plt.figure()
plt.plot(time[:datalen],farms_slip_z,'g-',label=farms_label)
plt.plot(benchmark[:,0],benchmark[:,4],'r-',label=benchmark_label)
plt.title(benchmark_code+" slip time history at strike "+str(xcoord_i)+"km and at dip "+str(zcoord_i)+"km ")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("slip (m)")
# plt.savefig(saved_path+"/slip.png")
plt.show()

## traction z
plt.figure()
plt.plot(time[:datalen],farms_traction_z/1e6,'g-',label=farms_label)
plt.plot(benchmark[:,0],benchmark[:,-1],'r-',label=benchmark_label)
plt.title(benchmark_code+" slip time history at strike "+str(xcoord_i)+"km and at dip "+str(zcoord_i)+"km ")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("slip (m)")
# plt.savefig(saved_path+"/slip.png")
plt.show()