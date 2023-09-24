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
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_920_dampingcomp/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_921_damp05/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/50m_921_damp05/"
# overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/100m_922_uguca_setup/"
overall_file_path = "/Volumes/One Touch/Research/RateStateDebug/2D/50m_922_uguca_setup/"

#exodus_file_path = overall_file_path + "main_out.e"
exodus_file_path = overall_file_path + "main_out.e"

# save_folder_output_file_path = "./files/100m"
# save_folder_output_file_path = "./files/50m"
# save_folder_output_file_path = "./files/100m_919"
# save_folder_output_file_path = "./files/100m_919_2"
# save_folder_output_file_path = "./files/50m_919"
# save_folder_output_file_path = "./files/100m_919_damp03"
# save_folder_output_file_path = "./files/100m_919_damp08"
# save_folder_output_file_path = "./files/100m_919_fixnormal"
# save_folder_output_file_path = "./files/50m_919_fixnormal"
# save_folder_output_file_path = "./files/100m_921_damp05"
# save_folder_output_file_path = "./files/100m_922_uguca_setup"
save_folder_output_file_path = "./files/50m_922_uguca_setup"

decodeflag = "name_nod_var"

#mid ptr locs index
#100m
# ptr1_loc = 51  - 1 #-5050
# ptr2_loc = 102 - 1 #50
# ptr3_loc = 152 - 1 #5050
# ptr4_loc = 76  - 1 #-2550
# ptr5_loc = 127 - 1 #2550
# ptr6_loc = 26  - 1 #-7550
# ptr7_loc = 177 - 1 #7550

#50m
ptr1_loc = 101 - 1 #-5025(101) #4975(102) #5075(103)
ptr2_loc = 202 - 1 #25
ptr3_loc = 302 - 1 #5025
ptr4_loc = 151 - 1 #-2525
ptr5_loc = 252 - 1 #2525
ptr6_loc = 51  - 1 #-7525
ptr7_loc = 352 - 1 #7525

#dt
dt = 0.025
t_max = 5.0
time = np.arange(0,t_max+dt,dt)
#save
np.savetxt(save_folder_output_file_path + "/time.txt",time)

#options
decodenow  = False
getvelnow = False
getslipratenow = False
plotnow = False

getdispnow = False
getresidnow = False
peakslipratepositionfunc = False
rupturespeed = True

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
    np.savetxt(save_folder_output_file_path + "/sliprate_neg5kmstrike0dot75dip.txt",sliprate[:,ptr1_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_0strike0dot75dip.txt",sliprate[:,ptr2_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_pos5kmstrike0dot75dip.txt",sliprate[:,ptr3_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_neg2dot5kmstrike0dot75dip.txt",sliprate[:,ptr4_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_pos2dot5kmstrike0dot75dip.txt",sliprate[:,ptr5_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_neg7dot5kmstrike0dot75dip.txt",sliprate[:,ptr6_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/sliprate_pos7dot5kmstrike0dot75dip.txt",sliprate[:,ptr7_loc])

if getdispnow == True:

    #get dispx
    dispx = nc.variables["vals_nod_var3"]

    dispx_upper = dispx[:,upper_ptr_index]
    dispx_lower = dispx[:,lower_ptr_index]

    slip = np.subtract(dispx_upper,dispx_lower)

    #save
    np.savetxt(save_folder_output_file_path + "/slip.txt",slip)
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_neg5kmstrike0dot75dip.txt",slip[:,ptr1_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_0strike0dot75dip.txt",slip[:,ptr2_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_pos5kmstrike0dot75dip.txt",slip[:,ptr3_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_neg2dot5kmstrike0dot75dip.txt",slip[:,ptr4_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_pos2dot5kmstrike0dot75dip.txt",slip[:,ptr5_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_neg7dot5kmstrike0dot75dip.txt",slip[:,ptr6_loc])
    #save middle ptrs time history
    np.savetxt(save_folder_output_file_path + "/slip_pos7dot5kmstrike0dot75dip.txt",slip[:,ptr7_loc])

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

if peakslipratepositionfunc == True:

    #get slip rate time series
    sliprate = np.loadtxt(save_folder_output_file_path + "/sliprate.txt")

    #load location
    interface_coordx = np.loadtxt(save_folder_output_file_path + "/interface_coordx.txt")

    sliprate_max = np.max(sliprate[:161,:],axis=0)

    #get vel0
    vel0 = np.loadtxt("/Users/andyz/projects/farms/postprocess/2d/TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx1440_s2.00_tf0.35_npc1-DataFiles/vel0.txt")
    sliprate0 = vel0 * 2
    sliprate_max0 = np.max(sliprate0[:40,:],axis=0)
    xcoord0 = np.linspace(-35975,35975,720*2)

    plt.plot(interface_coordx, sliprate_max, 'b.',label="moose-50m")
    plt.plot(xcoord0, sliprate_max0, 'r.',label="uguca-50m")
    plt.xlim([-10000,10000])
    plt.xlabel("x coordinate (m)")
    plt.ylabel("peak slip rate (m/s)")
    plt.legend()
    plt.show()

if rupturespeed == True:

    #load time
    time = np.loadtxt(save_folder_output_file_path + "/time.txt")

    #get slip rate time series
    sliprate = np.loadtxt(save_folder_output_file_path + "/sliprate.txt")

    #load location
    interface_coordx = np.loadtxt(save_folder_output_file_path + "/interface_coordx.txt")

    #get time
    timeright = time[:161]

    #get slip rate right half
    sliprateright = sliprate[:161,int(np.shape(sliprate)[1]/2):]

    #get location right half
    interface_coordxright = interface_coordx[int(np.shape(sliprate)[1]/2):]

    #sliprate threshold
    threshold = 0.1

    locs_plot = []
    time_plot = []
    index_save = []

    for i in range(np.shape(sliprateright)[0]):

        sliprateright_timei = sliprateright[i,:]
        
        index_rightmost = np.where(sliprateright_timei > threshold)[0].tolist()

        if np.shape(index_rightmost)[0] >= 1:
            index_rightmost = index_rightmost[-1]
        else:
            index_rightmost = 0

        if i != 0:
            if interface_coordxright[index_rightmost] < locs_plot[-1]:
                index_rightmost = index_save[-1]
        else:
            index_rightmost = 0

        locs_plot.append(interface_coordxright[index_rightmost])

        time_plot.append(timeright[i])

        index_save.append(index_rightmost)

    #fit whole curve
    polydata001_terms = np.polyfit(time_plot,locs_plot,6)
    polydata001 = np.poly1d(polydata001_terms)
    polyvel001_terms = np.polyder(polydata001_terms)
    polyvel001 = np.poly1d(polyvel001_terms)

    # plt.plot(time_plot,polydata001(time_plot))
    # plt.plot(time_plot,locs_plot)
    # plt.title("list_front_time0 vs list_front_xcoord0 or polydata001(list_front_time0)")
    # plt.show()

    # plt.plot(time_plot,polyvel001(time_plot))
    # plt.title("list_front_time0 vs polyvel001(list_front_time0)")
    # plt.show()

    ####for uguca data
    time = np.loadtxt("/Users/andyz/projects/farms/postprocess/2d/TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx1440_s2.00_tf0.35_npc1.time")
    time = time[:,1]
    #get slip rate time series
    vel = np.loadtxt("/Users/andyz/projects/farms/postprocess/2d/TPV101_Nx/TPV101_Nx_50and25m/TPV101_Nx1440_s2.00_tf0.35_npc1-DataFiles/vel0.txt")
    sliprate = vel * 2
    #load location
    interface_coordx = np.linspace(-35975,35975,720*2)

    #get time
    timeright = time[:40]
    #get slip rate right half
    sliprateright = sliprate[:40,int(np.shape(sliprate)[1]/2):]
    #get location right half
    interface_coordxright = interface_coordx[int(np.shape(sliprate)[1]/2):]

    #sliprate threshold
    threshold = 0.1
    ugucalocs_plot = []
    ugucatime_plot = []
    ugucaindex_save = []
    for i in range(np.shape(sliprateright)[0]):
        sliprateright_timei = sliprateright[i,:]
        index_rightmost = np.where(sliprateright_timei > threshold)[0].tolist()
        if np.shape(index_rightmost)[0] >= 1:
            index_rightmost = index_rightmost[-1]
        else:
            index_rightmost = 0

        if i != 0:
            if interface_coordxright[index_rightmost] < ugucalocs_plot[-1]:
                index_rightmost = ugucaindex_save[-1]
        else:
            index_rightmost = 0

        ugucalocs_plot.append(interface_coordxright[index_rightmost])
        ugucatime_plot.append(timeright[i])
        ugucaindex_save.append(index_rightmost)

    #fit whole curve
    ugucadata001_terms = np.polyfit(ugucatime_plot,ugucalocs_plot,6)
    ugucadata001 = np.poly1d(ugucadata001_terms)
    ugucavel001_terms = np.polyder(ugucadata001_terms)
    ugucavel001 = np.poly1d(ugucavel001_terms)

    # plt.plot(ugucatime_plot,ugucadata001(ugucatime_plot),'b-')
    plt.plot(ugucatime_plot,ugucalocs_plot,'r-.',label="uguca-50m")
    # plt.plot(time_plot,polydata001(time_plot),'r-')
    plt.plot(time_plot,locs_plot,'b-.',label="moose-50m")
    plt.ylim([0,10000])
    plt.xlim([0,4])
    plt.title("rupture tip position (m) vs time (s)")
    plt.xlabel("time (s)")
    plt.ylabel("rupture tip position (m)")
    plt.legend()
    plt.show()

    plt.plot(ugucatime_plot,ugucavel001(ugucatime_plot),'b--')
    plt.plot(time_plot,polyvel001(time_plot),'r--')
    plt.title("list_front_time0 vs polyvel001(list_front_time0)")
    plt.show()



if plotnow == True:

    #plot
    prop = np.loadtxt(save_folder_output_file_path + "/sliprate.txt")

    print(np.shape(prop))
    plt.figure()
    plt.plot(interface_coordx,prop[100-3,:])
    plt.show()




