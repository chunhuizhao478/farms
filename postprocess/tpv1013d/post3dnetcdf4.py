import numpy as np
import matplotlib.pyplot as plt

##loop time
num_begin = 0
num_end = 30
file_nums = np.linspace(num_begin,num_end,num_end+1)

time_begin = 0
time_end = 1.5
time = np.linspace(time_begin,time_end,num_end+1)

#data loc
loc = 502 - 1

#read data
# ##benchmark
# data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_DFM.txt")
# data_index_end = int(750) #DFM

# data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_liu.txt")
# data_index_end = int(156) #liu

# data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_MDSBI.txt")
# data_index_end = int(301) #MDSBI

# data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_SPECFEM3D.txt")
# data_index_end = int(300) #SPECFEM3D

# data_FaultMod100m = np.loadtxt("./benchmark/tpv102_0strike7dot5dip_Pylith_200m.txt")
# data_index_end = int(31) #Pylith

# data_FaultMod100m = np.loadtxt("./benchmark/tpv102_0strike7dot5dip_Pylith_200m_Tet4.txt")
# data_index_end = int(31) #Pylith

data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_uguca.txt")
data_index_end = int(149) #uguca

# data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_FaultMod.txt")
# data_index_end = int(150) #FaultMod

timehist_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,0]
traction_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,3]
sliprate_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,2]
statevar_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,-1]
slip_strikeFaultMod100m     = data_FaultMod100m[:data_index_end+1,1]

#load file
time_1 = np.loadtxt("./files/200m_fixnormal/time.txt")
sliprate_time_1 = np.loadtxt("./files/200m_fixnormal/sliprate_0strike0dot75dip.txt")
slip_time_1     = np.loadtxt("./files/200m_fixnormal/slip_0strike0dot75dip.txt")

time_2 = np.loadtxt("./files/100m_fixnormal/time.txt")
sliprate_time_2 = np.loadtxt("./files/100m_fixnormal/sliprate_0strike0dot75dip.txt")
slip_time_2     = np.loadtxt("./files/100m_fixnormal/slip_0strike0dot75dip.txt")
# statevar_time_2 = np.loadtxt("./files/100m_918/statevar_log_0strike0dot75dip.txt")
# traction_time_2 = np.loadtxt("./files/100m_918/traction_strike_0strike0dot75dip.txt")

# time_3          = np.loadtxt("./files/100m_testSep17/time.txt")
# sliprate_time_3 = np.loadtxt("./files/200m_918/sliprate_0strike0dot75dip.txt")
# slip_time_3     = np.loadtxt("./files/200m_918/slip_0strike0dot75dip.txt")
# statevar_time_3 = np.loadtxt("./files/200m_918/statevar_log_0strike0dot75dip.txt")
# traction_time_3 = np.loadtxt("./files/200m_918/traction_strike_0strike0dot75dip.txt")

#save data
plt.figure()
plt.plot(time_1,sliprate_time_1,'m-.',label="FARMS3D-200m-ApplyBC")
plt.plot(time_2,sliprate_time_2,'b-.',label="FARMS3D-100m-ApplyBC")
# plt.plot(time_3,sliprate_time_3,'k-.',label="FARMS3D-200m-ApplyBC")
plt.plot(timehist_strikeFaultMod100m,sliprate_strikeFaultMod100m,'r.',label='uguca-100m')
plt.title("Time History Slip Rate Along Strike Direction (m/s)")
plt.xlabel("time (s)")
plt.ylabel("slip rate (m/s)")
plt.legend()

#save data
plt.figure()
plt.plot(time_1,slip_time_1,'m-.',label="FARMS3D-200m-ApplyBC")
plt.plot(time_2,slip_time_2,'b-.',label="FARMS3D-100m-ApplyBC")
# plt.plot(time_3,slip_time_3,'k-.',label="FARMS3D-200m-ApplyBC")
plt.plot(timehist_strikeFaultMod100m,slip_strikeFaultMod100m,'r.',label='uguca-100m')
plt.title("Time History Slip Along Strike Direction (m/s)")
plt.xlabel("time (s)")
plt.ylabel("slip (m)")
plt.legend()

# #save data
# plt.figure()
# # plt.plot(time_1,sliprate_time_1,'m-.',label="FARMS3D-200m-ApplyBC")
# plt.plot(time_2,statevar_time_2,'b-.',label="FARMS3D-100m-ApplyBC")
# plt.plot(time_3,statevar_time_3,'k-.',label="FARMS3D-200m-ApplyBC")
# plt.plot(timehist_strikeFaultMod100m,statevar_strikeFaultMod100m,'r.',label='uguca-100m')
# plt.title("Time History State Variable Along Strike Direction (m/s)")
# plt.xlabel("time (s)")
# plt.ylabel("slip (m)")
# plt.legend()

# #save data
# plt.figure()
# # plt.plot(time_1,sliprate_time_1,'m-.',label="FARMS3D-200m-ApplyBC")
# plt.plot(time_2,traction_time_2,'b-.',label="FARMS3D-100m-ApplyBC")
# plt.plot(time_3,traction_time_3,'k-.',label="FARMS3D-200m-ApplyBC")
# plt.plot(timehist_strikeFaultMod100m,traction_strikeFaultMod100m,'r.',label='uguca-100m')
# plt.title("Time History Traction Along Strike Direction (m/s)")
# plt.xlabel("time (s)")
# plt.ylabel("slip (m)")
# plt.legend()

plt.show()