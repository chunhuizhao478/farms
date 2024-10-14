import numpy as np
import matplotlib.pyplot as plt

#coordinates
xcoord_i = 0.0
zcoord_i = 0.0

#time farms 0.1s interval
time = np.linspace(0,12.0,121)

#benchmark code
benchmark_code = "TPV205"

#read benchmark data
benchmark_label = "benchmark-eqsim-200m"

#read farms data
farms_label = "farms-200m"

#file path
benchmark_path = "./benchmark_data/eqsim_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
farms_path_Ty = "./farms_elem_data/Ty_strike0.0_dip0.0.txt"

benchmark = np.loadtxt(benchmark_path)
farms_Ty = np.loadtxt(farms_path_Ty,skiprows=1)

## check size of data
datalen = len(farms_Ty)

## slip rate
plt.figure()
plt.plot(time[:datalen],farms_Ty/1e6,'g-',label=farms_label)
plt.plot(benchmark[:,0],benchmark[:,6],'r-',label=benchmark_label)
plt.title(benchmark_code+" Traction y time history at strike "+str(xcoord_i)+"km and at dip "+str(zcoord_i)+"km ")
plt.legend()
plt.xlabel("time (s)")
plt.ylabel("slip rate (m/s)")
# plt.savefig(saved_path+"/sliprate.png")
plt.show() 
