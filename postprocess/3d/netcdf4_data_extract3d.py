import netCDF4
import numpy as np
import functions
import matplotlib.pyplot as plt

# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/3D/BC/200m_dipleak/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/3D/BC/100m_dipleak/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/3D/BC/100m_testSep17/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/3D/BC/50m/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/3D/BC/200m_918/"
overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/3D/BC/100m_918/"
exodus_file_path = overall_file_path + "main_out.e"

# save_folder_output_file_path = "./files/200m"
# save_folder_output_file_path = "./files/100m"
# save_folder_output_file_path = "./files/100m_testSep17"
# save_folder_output_file_path = "./files/50m"
# save_folder_output_file_path = "./files/200m_918"
save_folder_output_file_path = "./files/100m_918"

decodeflag = "name_nod_var"

#mid ptr locs
mid_ptr_loc = 50

#dt
dt = 0.05
t_max = 1.5
num_end = 30
time = np.arange(0,t_max+dt,dt)

#file nums
num_begin = 0
num_end = 30
file_nums = np.linspace(num_begin,num_end,num_end+1)

#save
np.savetxt(save_folder_output_file_path + "/time.txt",time)

#options
decodenow  = False
getvelnow = False
getslipratenow = False
plotnow = True

getdispnow = False
getslipnow = False

get_traction_statevar_now = False

nc = netCDF4.Dataset(exodus_file_path)

#!check all variables names
print(nc.variables.keys())

coordx = nc.variables["coordx"]
coordy = nc.variables["coordy"]

#read points data
upper_ptr_data = np.loadtxt(overall_file_path+"block0_block1_ptrs.txt",skiprows=1,delimiter=',')
lower_ptr_data = np.loadtxt(overall_file_path+"block1_block0_ptrs.txt",skiprows=1,delimiter=',')

#get specific dip location points
upper_ptr_data_zcoord = upper_ptr_data[:,2]
upper_ptr_data_dip0_index = np.where((upper_ptr_data_zcoord > -1) & (upper_ptr_data_zcoord < 1))[0]
# upper_ptr_data_dip0_index = np.where((upper_ptr_data_zcoord > -7550) & (upper_ptr_data_zcoord < -7450))[0]
upper_ptr_data = upper_ptr_data[upper_ptr_data_dip0_index,:]

lower_ptr_data_zcoord = lower_ptr_data[:,2]
lower_ptr_data_dip0_index = np.where((lower_ptr_data_zcoord > -1) & (lower_ptr_data_zcoord < 1))[0]
# lower_ptr_data_dip0_index = np.where((lower_ptr_data_zcoord > -7550) & (lower_ptr_data_zcoord < -7450))[0]
lower_ptr_data = lower_ptr_data[lower_ptr_data_dip0_index,:]

#get x coordinate
interface_coordx = upper_ptr_data[:,0]

#print z coordinate
print("z coordinate: ", upper_ptr_data[:,2][0])

#get index (ids - 1)
upper_ptr_index = upper_ptr_data[:,-1] - 1
lower_ptr_index = lower_ptr_data[:,-1] - 1

#convert into integer
upper_ptr_index = [int(i) for i in upper_ptr_index]
lower_ptr_index = [int(i) for i in lower_ptr_index]

#!check all points on the interface
# for i in lower_ptr_index:
#     y_i = coordy[i]
#     if y_i != 0:
#         print("wrong")

traction_strike = []
statevar_log = []

#Decode Nodal Name
if decodenow == True:   
    
    functions.DecodeNameNodal(nc, decodeflag, save_folder_output_file_path)

if getvelnow == True:

    #get vel_x
    velx = nc.variables["vals_nod_var10"]
    #save velx
    np.savetxt(save_folder_output_file_path + "/velx.txt",velx)

if getslipratenow == True:

    velx = np.loadtxt(save_folder_output_file_path + "/velx.txt")

    upper_velx = velx[:,upper_ptr_index]
    lower_velx = velx[:,lower_ptr_index]

    #elementwise ops
    sliprate = np.subtract(upper_velx,lower_velx)
    #save
    np.savetxt(save_folder_output_file_path + "/sliprate.txt",sliprate)
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_0strike0dot75dip.txt",sliprate[:,mid_ptr_loc])

if getdispnow == True:

    #get dispx
    dispx = nc.variables["vals_nod_var4"]
    #save dispx
    np.savetxt(save_folder_output_file_path + "/dispx.txt",dispx)

if getslipnow == True:

    dispx = np.loadtxt(save_folder_output_file_path + "/dispx.txt")

    upper_dispx = dispx[:,upper_ptr_index]
    lower_dispx = dispx[:,lower_ptr_index]

    #elementwise ops
    slip = np.subtract(upper_dispx,lower_dispx)
    #save
    np.savetxt(save_folder_output_file_path + "/slip.txt",slip)
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_0strike0dot75dip.txt",slip[:,mid_ptr_loc])

if get_traction_statevar_now == True:

   loc = 502 - 1

   for index_i in range(num_end+1):
       
       #get file 
       file_num_i = file_nums[index_i]
       
       #load file
       datafile_time = np.loadtxt(save_folder_output_file_path+"/datatimeseries/data_3d_100m_"+str(int(file_num_i))+".txt", skiprows=1, delimiter=',')
       
       #get dataline
       data_time = datafile_time[loc,:]

       #save data
       traction_strike.append(data_time[2]/1e6)
       statevar_log.append(np.log10(data_time[1]))
    
   #save file
   np.savetxt(save_folder_output_file_path + "/traction_strike_0strike0dot75dip.txt",traction_strike)
   np.savetxt(save_folder_output_file_path + "/statevar_log_0strike0dot75dip.txt",statevar_log)

if plotnow == True:

    #plot
    sliprate = np.loadtxt(save_folder_output_file_path + "/sliprate.txt")

    for i in range(num_end):
    # print(np.shape(interface_coordx))
        plt.figure()
        plt.plot(interface_coordx,sliprate[i,:])
        plt.ylim([0,4.5])
        plt.savefig(save_folder_output_file_path+'/sliprate'+str(i)+'.png')
        plt.close()




