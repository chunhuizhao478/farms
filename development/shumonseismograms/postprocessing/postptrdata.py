import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.spatial import distance

def DecodeName(nc, decodeflag, save_folder_output_file_path):

    """
    Get Nodal Var Name 
    """

    # #Read File
    # nc = netCDF4.Dataset(exodus_file_path)

    #Get numpy.bytes name array
    if decodeflag == "name_nod_var":
        arr_name_nod_var = nc.variables['name_nod_var']

    #Get number of name
    num_names = np.shape(arr_name_nod_var)[0]

    #initialize list of names
    list_name_nod_var = []

    #loop over each name
    for name_ind in range(num_names):

        name_i_decode = ''

        name_i = arr_name_nod_var[name_ind]

        # print(name_i)

        for ind in range(len(name_i)):
            
            if not np.ma.is_masked(name_i[ind]):
                
                byte_i_decode = name_i[ind].decode('UTF-8')

                name_i_decode += byte_i_decode

        list_name_nod_var.append(name_i_decode)
    
    #save
    if decodeflag == "name_nod_var":
        np.savetxt(save_folder_output_file_path + "/list_name_nod_var.txt",list_name_nod_var,fmt='%s',newline=" ")

def MaxVal(nc, 
           save_folder_output_file_path, 
           plot_var_name):

    ##
    #Initialize dict storing elem var name / index pair
    dict_evn_index = {}

    #Read Elem Var Name
    evn = np.loadtxt(save_folder_output_file_path + "/list_name_nod_var.txt",dtype=str)

    #fill in dict_evn_index
    for evn_ind in range(len(evn)):
        evn_i = evn[evn_ind]
        # start from 1
        dict_evn_index[evn_i] = evn_ind + 1

    #get required param to be plotted
    #loop over plot names
    var_ind_list = []
    for plot_var_name_index in range(len(plot_var_name)):

        plot_var_name_i = plot_var_name[plot_var_name_index]
        var_ind = dict_evn_index[plot_var_name_i]

        var_ind_list.append(var_ind)

    # print(dict_evn_index)
    # print(var_ind_list)
    # exit()

    #get data_x, data_y #hardcode 2D, only vel#
    #(num of time steps, num of nodes)
    arr_param_x = nc.variables['vals_nod_var'+str(var_ind_list[0])][:,:]
    arr_param_y = nc.variables['vals_nod_var'+str(var_ind_list[1])][:,:]

    #elementwise array operations
    arr_param_x_power2 = np.power(arr_param_x,2)
    arr_param_y_power2 = np.power(arr_param_y,2)
    arr_param_addtwo = np.add(arr_param_x_power2,arr_param_y_power2)
    arr_param_sqrt   = np.sqrt(arr_param_addtwo)

    # print(np.shape(arr_param_sqrt))
    # exit()
    
    #find maximum value across time at each node
    arr_param_tmax = np.amax(arr_param_sqrt,axis=0)

    # print(np.shape(arr_param_tmax))

    #load coordinate
    x_coord = nc.variables['coordx']
    y_coord = nc.variables['coordy']

    np.savetxt(save_folder_output_file_path + "/maxvelvals.txt",arr_param_tmax,fmt='%.3f',newline=" ")
    np.savetxt(save_folder_output_file_path + "/xcoordvals.txt",x_coord[:],fmt='%.3f',newline=" ")
    np.savetxt(save_folder_output_file_path + "/ycoordvals.txt",y_coord[:],fmt='%.3f',newline=" ")

def PlotPtrXYVal(nc,save_folder_output_file_path,plot_var_name,ptr_coord,angle,i,showfig=False):

    #plot x y component values for given "plot_var_name" at given point "ptr_coord"

    ##
    #Initialize dict storing elem var name / index pair
    dict_evn_index = {}

    #Read Elem Var Name
    evn = np.loadtxt(save_folder_output_file_path + "/list_name_nod_var.txt",dtype=str)

    #fill in dict_evn_index
    for evn_ind in range(len(evn)):
        evn_i = evn[evn_ind]
        # start from 1
        dict_evn_index[evn_i] = evn_ind + 1

    #get required param to be plotted
    #loop over plot names
    var_ind_list = []
    for plot_var_name_index in range(len(plot_var_name)):

        plot_var_name_i = plot_var_name[plot_var_name_index]
        var_ind = dict_evn_index[plot_var_name_i]

        var_ind_list.append(var_ind)
    
    print("loading val array ...")

    #get data_x, data_y #hardcode 2D, only vel#
    #(num of time steps, num of nodes)
    arr_param_x = nc.variables['vals_nod_var'+str(var_ind_list[0])][:,:]
    arr_param_y = nc.variables['vals_nod_var'+str(var_ind_list[1])][:,:]

    #current point x y coordinate
    ptr_x = ptr_coord[0]
    ptr_y = ptr_coord[1]

    print("loading coordinate array ...")

    #load coordinate
    x_coord = nc.variables['coordx']
    y_coord = nc.variables['coordy']

    #find index of that point
    idc = helper_find_nearest(x_coord[:], ptr_x, y_coord[:], ptr_y)

    # print(idc)
    # print(np.shape(idc))
    # exit()
    #get time

    print("loading time array ...")

    #Obtain Time Series
    timeseries = nc.variables['time_whole']

    print("slicing the x,y values ...")

    #get time history of current ptr
    ptr_x_valhist = arr_param_x[:,idc]
    ptr_y_valhist = arr_param_y[:,idc]

    # print("compute the magnitude of given variable values ...")

    # ptr_mag_valhist = np.sqrt(ptr_x_valhist**2+ptr_y_valhist**2)

    # print("direction converting to fault normal/parallel ...")
    
    # #fault_normal/fault_parallel components
    # theta = math.radians(angle)
    # #fault_normal = np.multiply(math.sin(theta),ptr_x_valhist) + np.multiply(math.cos(theta),ptr_y_valhist)
    # #fault_parallel = np.multiply(math.cos(theta),ptr_x_valhist) - np.multiply(math.sin(theta),ptr_y_valhist)

    # fault_normal = np.multiply(math.sin(theta),ptr_x_valhist) - np.multiply(math.cos(theta),ptr_y_valhist)
    # fault_parallel = np.multiply(math.cos(theta),ptr_x_valhist) + np.multiply(math.sin(theta),ptr_y_valhist)

    # print("filtering . remove high frenquency ...")

    # create a normalized Hanning window
    # windowSize = 40
    # window = np.hanning(windowSize)
    # window = window / window.sum()

    # # filter the data using convolution
    # faultnormalfiltered = np.convolve(window, fault_normal, mode='same')
    # faultparallelfiltered = np.convolve(window, fault_parallel, mode='same')

    print("Performing Fourier Transform ...") 

    # Perform the Fourier Transform
    ## time
    frequencies = np.fft.fftfreq(timeseries[:].size, timeseries[:][1] - timeseries[:][0])   

    ## x dir
    V_fft_x = np.fft.fft(ptr_x_valhist)
    
    # We need only the positive half of the spectrum, since the negative is a mirror of the positive
    pos_half = frequencies > 0
    
    # Find the peak magnitude and its corresponding frequency
    peak_magnitude_x = np.max(np.abs(V_fft_x[pos_half]))
    peak_frequency_x = frequencies[pos_half][np.argmax(np.abs(V_fft_x[pos_half]))]
    
    ## y dir
    V_fft_y = np.fft.fft(ptr_y_valhist)

    # We need only the positive half of the spectrum, since the negative is a mirror of the positive
    pos_half = frequencies > 0
    
    # Find the peak magnitude and its corresponding frequency
    peak_magnitude_y = np.max(np.abs(V_fft_y[pos_half]))
    peak_frequency_y = frequencies[pos_half][np.argmax(np.abs(V_fft_y[pos_half]))]

    print("ploting ...")

    #plot velocity in x
    #global
    plt.figure(figsize=(10,6))
    plt.subplot(1, 2, 1)
    plt.plot(timeseries[:],ptr_x_valhist,'b-',label = "velocity in x direction")
    # plt.plot(timeseries[:],ptr_y_valhist,'r-',label = "velocity y")
    plt.xlabel("time history (s)")
    plt.ylabel("velocity (m/s)")
    plt.legend()
    plt.title("Time History of Velocity at location: "+str(ptr_x)+" , "+str(ptr_y))
    # plt.xlim([0,35])
    # plt.ylim([0,18.0])

    plt.subplot(1, 2, 2)  
    plt.plot(frequencies, np.abs(V_fft_x))
    plt.title('Magnitude of Fourier Transform')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude')     
    plt.tight_layout()

    plt.savefig(save_folder_output_file_path+'/imgx'+str(i)+'.png')

    #plot velocity in y
    #global
    plt.figure(figsize=(10,6))
    plt.subplot(1, 2, 1)
    plt.plot(timeseries[:],ptr_y_valhist,'b-',label = "velocity in y direction")
    # plt.plot(timeseries[:],ptr_y_valhist,'r-',label = "velocity y")
    plt.xlabel("time history (s)")
    plt.ylabel("velocity (m/s)")
    plt.legend()
    plt.title("Time History of Velocity at location: "+str(ptr_x)+" , "+str(ptr_y))
    # plt.xlim([0,35])
    # plt.ylim([0,18.0]) 
    
    plt.subplot(1, 2, 2)  
    plt.plot(frequencies, np.abs(V_fft_y))
    plt.title('Magnitude of Fourier Transform')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude') 
    plt.tight_layout()

    plt.savefig(save_folder_output_file_path+'/imgy'+str(i)+'.png')
    
    #
    np.savetxt(save_folder_output_file_path+'/list_timeseries'+str(i)+'.txt',
                    timeseries,
                    fmt='%.7f',
                    newline=" ")
    np.savetxt(save_folder_output_file_path+'/list_velx'+str(i)+'.txt',
                    ptr_x_valhist,
                    fmt='%.7f',
                    newline=" ")
    np.savetxt(save_folder_output_file_path+'/list_vely'+str(i)+'.txt',
                    ptr_y_valhist,
                    fmt='%.7f',
                    newline=" ")
    np.savetxt(save_folder_output_file_path+'/list_peaks'+str(i)+'.txt',
                    [peak_magnitude_x,peak_frequency_x,peak_magnitude_y,peak_frequency_y],
                    fmt='%.7f',
                    newline=" ")
    

    # #fault_local
    # plt.figure()
    # plt.plot(timeseries[:],fault_normal,'b-',label = "fault normal")
    # plt.plot(timeseries[:],fault_parallel,'r-',label = "fault parallel")
    # plt.xlabel("time history (s)")
    # plt.ylabel("velocity (m/s)")
    # plt.legend()
    # plt.title("Time History of Velocity (Fault)")
    # #splay fault
    # # plt.xlim([0,15])

    # #fault_local_filtered
    # plt.figure()
    # plt.plot(timeseries[:],faultnormalfiltered,'b-',label = "fault normal")
    # plt.plot(timeseries[:],faultparallelfiltered,'r-',label = "fault parallel")
    # plt.xlabel("time history (s)")
    # plt.ylabel("velocity (m/s)")
    # plt.legend()
    # plt.title("Time History of Velocity (Fault Filtered)")

    if showfig:
        plt.show()

def helper_find_nearest(array_x, value_x, array_y, value_y):
    
    print("start finding nearest point ...")

    #
    array_x = array_x[:].reshape((len(array_x),1))
    array_y = array_y[:].reshape((len(array_y),1))

    #
    nodes = np.hstack((array_x,array_y))
    node = np.hstack((value_x,value_y))

    closest_index = distance.cdist([node], nodes).argmin()

    print("find point! x: "+str(nodes[closest_index][0])+" y: "+str(nodes[closest_index][1]))

    return closest_index

    # #x
    # array_x = np.asarray(array_x)
    # idx = (np.abs(array_x - value_x)).argmin()

    # print("finish finding nearest point index x coord ...")

    # #y
    # array_y = np.asarray(array_y)
    # idy = (np.abs(array_y - value_y)).argmin()

    # print("finish finding nearest point index y coord ...")

    #find common values
    # idc = np.intersect1d(idx,idy)[0] 

    # if abs(value_x - array_x[idx]) < 1 and abs(value_y - array_y[idx]) < 1:
    #     print("find point! x: "+str(array_x[idx])+" y: "+str(array_y[idx]))
    #     return idx
    # elif abs(value_x - array_x[idy]) < 1 and abs(value_y - array_y[idy]) < 1:
    #     print("find point! y: "+str(array_x[idy])+" y: "+str(array_y[idy]))
    #     return idy
    # else:
    #     print("unable to find point! Exit ...")
    #     exit()
    #     return None


def PlotMaxVal(nc,save_folder_output_file_path):

    #load val
    arr_param = np.loadtxt(save_folder_output_file_path + "/maxvelvals.txt")

    #load coordinate
    x_coord = nc.variables['coordx']
    y_coord = nc.variables['coordy']

    #plot
    ax = plt.axes(projection='3d')
    ax.scatter3D(x_coord, y_coord, arr_param, c=arr_param, cmap='PuBu');
    plt.show()

def WriteNCFile(nc,save_folder_output_file_path):

    #load val
    arr_param = np.loadtxt(save_folder_output_file_path + "/maxvelvals.txt")

    #load coordinate
    x_coord = nc.variables['coordx']
    y_coord = nc.variables['coordy']

    #create nc dataset
    fn = save_folder_output_file_path + "/velout.nc"
    ds = netCDF4.Dataset(fn,'w',format="NETCDF4")

    #create dimension
    time = ds.createDimension('time', None) 
    x = ds.createDimension('x',np.shape(x_coord)[0])
    y = ds.createDimension('y',np.shape(y_coord)[0])

    #create variables
    times = ds.createVariable('time', 'f4', ('time',))
    xs = ds.createVariable('x', 'f4', ('x',))
    ys = ds.createVariable('y', 'f4', ('y',))
    value = ds.createVariable('value', 'f4', ('time', 'x', 'y',))
    value.units = 'Unknown'

    xs[:] = x_coord[:]
    ys[:] = y_coord[:]
    
    value[0, :, :] = arr_param[:]

    ds.close()



#file path
exodus_file_path = "../inputfiles/explicitdynamic_1_out_damped.e"
save_folder_output_file_path = "./outputs"

#read exodus file
nc = netCDF4.Dataset(exodus_file_path)

#decode name
decodeflag = "name_nod_var"

#plot variable names
plot_var_name = ["vel_x","vel_y"]

#redo get val op
get_val_op = False

given_coord_list = [[0, -180],
                    [0, -360],
                    [0, -540]]

#angle of faults
angle = 0 #deg #not used

# print(nc.variables.keys()) #name_nod_var
# print(nc.variables["vals_nod_var10"])

##Decode name
DecodeName(nc, decodeflag ,save_folder_output_file_path)

##Output Maximum Vals
if get_val_op == True:
    MaxVal(nc, 
        save_folder_output_file_path, 
        plot_var_name)

for i in range(np.shape(given_coord_list)[0]):

    given_coord = given_coord_list[i]

    #Plot Time History of Point Quantity
    PlotPtrXYVal(nc,save_folder_output_file_path,plot_var_name,given_coord,angle,i,showfig=True)

# Save points locations
# np.savetxt("./outputs/given_coord_list.txt",given_coord_list)

##Plot
# PlotMaxVal(nc,save_folder_output_file_path)

##Write File
# WriteNCFile(nc,save_folder_output_file_path)