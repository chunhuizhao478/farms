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


def PlotPtrVal(nc,save_folder_output_file_path,plot_var_name,ptr_coord,showfig=False):

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

    #get data_x, data_y
    #(num of time steps, num of nodes)
    arr_param_vel  = nc.variables['vals_nod_var'+str(var_ind_list[0])][:,:]
    arr_param_disp = nc.variables['vals_nod_var'+str(var_ind_list[1])][:,:]

    #current point x y coordinate
    ptr_x = ptr_coord[0]
    ptr_y = ptr_coord[1]
    ptr_z = ptr_coord[2]

    print("loading coordinate array ...")

    #load coordinate
    x_coord = nc.variables['coordx']
    y_coord = nc.variables['coordy']
    z_coord = nc.variables['coordz']

    #find index of that point
    idc = helper_find_nearest(x_coord[:], ptr_x, y_coord[:], ptr_y, z_coord[:], ptr_z)

    #get time

    print("loading time array ...")

    #Obtain Time Series
    timeseries = nc.variables['time_whole']

    print("slicing the x,y values ...")

    #get time history of current ptr
    ptr_sliprate_valhist = 2 * abs( arr_param_vel[:,idc] )
    ptr_slip_valhist = 2 * abs( arr_param_disp[:,idc] )

    #
    np.savetxt(save_folder_output_file_path+'/sliprate_strike'+str(ptr_x/1000)+'_dip'+str(ptr_z/1000)+'.txt',
                ptr_sliprate_valhist,
                fmt='%.7f',
                newline=" ")
    np.savetxt(save_folder_output_file_path+'/slip_strike'+str(ptr_x/1000)+'_dip'+str(ptr_z/1000)+'.txt',
                ptr_slip_valhist,
                fmt='%.7f',
                newline=" ")

def helper_find_nearest(array_x, value_x, array_y, value_y, array_z, value_z):
    
    print("start finding nearest point ...")

    #
    array_x = array_x[:].reshape((len(array_x),1))
    array_y = array_y[:].reshape((len(array_y),1))
    array_z = array_z[:].reshape((len(array_z),1))

    #
    nodes = np.hstack((array_x,array_y,array_z))
    node = np.hstack((value_x,value_y,value_z))

    closest_index = distance.cdist([node], nodes).argmin()

    print("find point! x: "+str(nodes[closest_index][0])+" y: "+str(nodes[closest_index][1])+" z: "+str(nodes[closest_index][2]))

    return closest_index

#file path
exodus_file_path = "/Users/zhaoc/Downloads/tpv2053D_out.e"
save_folder_output_file_path = "./farms_data"

#read exodus file
nc = netCDF4.Dataset(exodus_file_path)

#decode name
decodeflag = "name_nod_var"

#plot variable names
plot_var_name = ["vel_slipweakening_x","disp_slipweakening_x"]

#strike,dip
given_coord_list = [[0     , 0 ,-0    ],
                    [0     , 0 ,-3000 ],
                    [0     , 0 ,-7500 ],
                    [0     , 0 ,-12000],
                    [4500  , 0 ,-0    ],
                    [4500  , 0 ,-7500 ],
                    [-4500 , 0 ,-0    ],
                    [-4500 , 0 ,-7500 ],
                    [7500  , 0 ,-0    ],
                    [7500  , 0 ,-7500 ],
                    [-7500 , 0 ,-0    ],
                    [-7500 , 0 ,-7500 ],
                    [12000 , 0 ,-0    ],
                    [12000 , 0 ,-7500 ],
                    [-12000, 0 ,-0    ],
                    [-12000, 0 ,-7500 ]]

##Decode name
DecodeName(nc, decodeflag ,save_folder_output_file_path)

for i in range(np.shape(given_coord_list)[0]):

    given_coord = given_coord_list[i]

    #Plot Time History of Point Quantity
    PlotPtrVal(nc,save_folder_output_file_path,plot_var_name,given_coord,showfig=False)

# Save points locations
np.savetxt("./outputs/given_coord_list.txt",given_coord_list)