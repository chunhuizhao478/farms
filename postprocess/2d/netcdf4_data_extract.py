import netCDF4
import numpy as np
import functions
import matplotlib.pyplot as plt

# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/50m/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_919/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_919_2/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/50m_919/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_919_damp03/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_919_damp08/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_919_fixnormal/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/50m_919_fixnormal/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/25m_919_fixnormal/"
overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_920_dampingcomp/"

#exodus_file_path = overall_file_path + "main_out.e"
exodus_file_path = overall_file_path + "main_out_1dot5damp.e"

# save_folder_output_file_path = "./files/100m"
# save_folder_output_file_path = "./files/50m"
# save_folder_output_file_path = "./files/100m_919"
# save_folder_output_file_path = "./files/100m_919_2"
# save_folder_output_file_path = "./files/50m_919"
# save_folder_output_file_path = "./files/100m_919_damp03"
# save_folder_output_file_path = "./files/100m_919_damp08"
# save_folder_output_file_path = "./files/100m_919_fixnormal"
# save_folder_output_file_path = "./files/50m_919_fixnormal"
save_folder_output_file_path = "./files/100m_920_dampcomp_1dot5"

decodeflag = "name_nod_var"

#mid ptr locs
mid_ptr_loc = 100

#dt
dt = 0.025
t_max = 3
time = np.arange(0,t_max+dt,dt)
#save
np.savetxt(save_folder_output_file_path + "/time.txt",time)

#options
decodenow  = True
getvelnow = True
getslipratenow = True
plotnow = True

getdispnow = True
getresidnow = True

nc = netCDF4.Dataset(exodus_file_path)

#!check all variables names
print(nc.variables.keys())

coordx = nc.variables["coordx"]
coordy = nc.variables["coordy"]

#read points data
upper_ptr_data = np.loadtxt(overall_file_path+"block0_block1_ptrs.txt",skiprows=1,delimiter=',')
lower_ptr_data = np.loadtxt(overall_file_path+"block1_block0_ptrs.txt",skiprows=1,delimiter=',')

#get x coordinate
interface_coordx = upper_ptr_data[:,0]

np.savetxt(save_folder_output_file_path + "/interface_coordx.txt",interface_coordx)

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

#Decode Nodal Name
if decodenow == True:   
    
    functions.DecodeNameNodal(nc, decodeflag, save_folder_output_file_path)

if getvelnow == True:

    #get vel_x
    velx = nc.variables["vals_nod_var9"]
    #save velx
    np.savetxt(save_folder_output_file_path + "/velx.txt",velx)

if getslipratenow == True:

    velx = np.loadtxt(save_folder_output_file_path + "/velx.txt")

    upper_velx = velx[:,upper_ptr_index]
    lower_velx = velx[:,lower_ptr_index]

    sliprate = np.zeros((np.shape(upper_velx)[0],np.shape(upper_velx)[1]))
    #elementwise ops
    sliprate = np.subtract(upper_velx,lower_velx)
    #save
    np.savetxt(save_folder_output_file_path + "/sliprate.txt",sliprate)
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_0strike0dot75dip.txt",sliprate[:,mid_ptr_loc])

if getdispnow == True:

    #get dispx
    dispx = nc.variables["vals_nod_var3"]

    dispx_upper = dispx[:,upper_ptr_index]
    np.savetxt(save_folder_output_file_path + "/dispx_upper.txt",dispx_upper)

    dispx_lower = dispx[:,upper_ptr_index]
    np.savetxt(save_folder_output_file_path + "/dispx_lower.txt",dispx_lower)

    #get dispy
    dispy = nc.variables["vals_nod_var4"]

    dispy_upper = dispy[:,upper_ptr_index]
    #save dispx
    np.savetxt(save_folder_output_file_path + "/dispy_upper.txt",dispy_upper)

if getresidnow == True:

    #get residx
    residx = nc.variables["vals_nod_var7"]

    residx_upper = residx[:,upper_ptr_index]
    #save residx
    np.savetxt(save_folder_output_file_path + "/residx_upper.txt",residx_upper)

if plotnow == True:

    #plot
    prop = np.loadtxt(save_folder_output_file_path + "/sliprate.txt")

    print(np.shape(prop))
    plt.figure()
    plt.plot(interface_coordx,prop[100-3,:])
    plt.show()




