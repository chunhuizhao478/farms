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

data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_uguca.txt")
data_index_end = int(150) #uguca

# data_FaultMod100m = np.loadtxt("./benchmark/tpv101_0strike7dot5dip_FaultMod.txt")
# data_index_end = int(150) #FaultMod

timehist_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,0]
traction_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,3]
sliprate_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,2]
statevar_strikeFaultMod100m = data_FaultMod100m[:data_index_end+1,-1]
slip_strikeFaultMod100m     = data_FaultMod100m[:data_index_end+1,1]

#initialize container
traction_strike_0 = []
sliprate_strike_0 = []
statevar_log_0 = []
slip_strike_0 = []

traction_strike = []
sliprate_strike = []
statevar_log = []
slip_strike = []

traction_strike_2 = []
sliprate_strike_2 = []
statevar_log_2 = []
slip_strike_2 = []

for index_i in range(num_end+1):

    #get file 
    file_num_i = file_nums[index_i]

    #load file
    datafile_time_h = np.loadtxt("./data/200m_dipleak/data_3d_200m_"+str(int(file_num_i))+".txt", skiprows=1, delimiter=',')
    datafile_time_i = np.loadtxt("./data/100m/data_3d_100m_"+str(int(file_num_i))+".txt", skiprows=1, delimiter=',')
    datafile_time_j = np.loadtxt("./data/50m/data_3d_50m_"+str(int(file_num_i))+".txt", skiprows=1, delimiter=',')

    #get dataline
    data_time_h = datafile_time_h[loc,:]
    data_time_i = datafile_time_i[loc,:]
    data_time_j = datafile_time_j[loc,:]

    #save data
    traction_strike_0.append(data_time_h[4]/1e6)
    statevar_log_0.append(np.log10(data_time_h[3]))
    sliprate_strike_0.append(data_time_h[2])
    slip_strike_0.append(data_time_h[1])

    traction_strike.append(data_time_i[4]/1e6)
    statevar_log.append(np.log10(data_time_i[3]))
    sliprate_strike.append(data_time_i[2])
    slip_strike.append(data_time_i[1])

    traction_strike_2.append(data_time_j[4]/1e6)
    statevar_log_2.append(np.log10(data_time_j[3]))
    sliprate_strike_2.append(data_time_j[2])
    slip_strike_2.append(data_time_j[1])

plt.figure()
plt.plot(time[1:],traction_strike_0[1:],'m-.',label="FARMS3D-200m-ApplyBC")
plt.plot(time[1:],traction_strike[1:],'b-.',label="FARMS3D-100m-ApplyBC")
plt.plot(time[1:],traction_strike_2[1:],'k-.',label="FARMS3D-50m-ApplyBC")
plt.plot(timehist_strikeFaultMod100m,traction_strikeFaultMod100m,'r.',label='uguca-100m')
plt.title("Time History Traction Along Strike Direction (MPa)")
plt.xlabel("time (s)")
plt.ylabel("traction (MPa)")
plt.legend()

plt.figure()
plt.plot(time[1:],sliprate_strike_0[1:],'m-.',label="FARMS3D-200m-ApplyBC")
plt.plot(time[1:],sliprate_strike[1:],'b-.',label="FARMS3D-100m-ApplyBC")
plt.plot(time[1:],sliprate_strike_2[1:],'k-.',label="FARMS3D-50m-ApplyBC")
plt.plot(timehist_strikeFaultMod100m,sliprate_strikeFaultMod100m,'r.',label='uguca-100m')
plt.title("Time History Slip Rate Along Strike Direction (m/s)")
plt.xlabel("time (s)")
plt.ylabel("slip rate (m/s)")
plt.legend()

plt.figure()
plt.plot(time[1:],statevar_log_0[1:],'m-.',label="FARMS3D-200m-ApplyBC")
plt.plot(time[1:],statevar_log[1:],'b-.',label="FARMS3D-100m-ApplyBC")
plt.plot(time[1:],statevar_log_2[1:],'k-.',label="FARMS3D-50m-ApplyBC")
plt.plot(timehist_strikeFaultMod100m,statevar_strikeFaultMod100m,'r.',label='uguca-100m')
plt.title("Time History log state variable")
plt.xlabel("time (s)")
plt.ylabel("state variable")
plt.legend()

plt.figure()
plt.plot(time[1:],slip_strike_0[1:],'m-.',label="FARMS3D-200m-ApplyBC")
plt.plot(time[1:],slip_strike[1:],'b-.',label="FARMS3D-100m-ApplyBC")
plt.plot(time[1:],slip_strike_2[1:],'k-.',label="FARMS3D-50m-ApplyBC")
plt.plot(timehist_strikeFaultMod100m,slip_strikeFaultMod100m,'r.',label='uguca-100m')
plt.title("Time History slip")
plt.xlabel("time (s)")
plt.ylabel("slip")
plt.legend()

#
plt.show()