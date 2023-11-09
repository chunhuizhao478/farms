"""
Post-Processing and Visualization
"""

#import module
import numpy as np
import meshio
import netCDF4
import matplotlib.pyplot as plt
import os 
import math
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable

"""
Files need to be generated in paraview:

- czm_ptrs_data.txt : file includes all points values along branch faults
0. de-select "Apply Displacements"
1. de-select all "Element Blocks", selects "Set" czm
Representation choose "Points"
2. File - Save Data - type name "czm_ptrs_data.txt"
3. Choose Arrays to Write: Array Selection:
GlobalNodeId, Point_x, Point_y, Point_z
Field Association: Point Data

Note: three files in current case:

1. mf180_ptrs_data.txt

"""

def RetrieveMFElementIDs(exodus_file_path):

    """
    Retrieve all Element IDs for Main Fault Region
    """

    ###--------------------------Read exodus file-----------------------------------###
    #read file
    nc = netCDF4.Dataset(exodus_file_path)

    #get number of element block
    num_elem_block = len(nc.variables['eb_status'])
    ###-----------------------------------------------------------------------------###

    #loop over element block
    for elemblock_ind in range(1,num_elem_block+1):

        #track progress
        print("(Main) elemblock_ind: ", elemblock_ind)
        print("(Main) Progress: ", elemblock_ind/(num_elem_block)*100,"%")

        #initialize list for saving elem along main fault
        list_mf_elem = []

        #initialize list for saving x coordinate
        list_mf_xcoord = []

        #initialize list for saving y coordinate
        list_mf_ycoord = []

        #initialize list for saving centroid point x coord of tria element (for ploting)
        list_mf_xcoordcenter = []

        #initialize list for saving centroid point y coord of tria element (for ploting)
        list_mf_ycoordcenter = []

        #get connectivity for the current block
        #shrift by 1
        tria_elem_connect = nc.variables['connect' + str(elemblock_ind)][:]-1

        #TRIA3 Elem Num
        num_elem = np.shape(tria_elem_connect)[0]

        #loop over element
        for elem_ind in range(num_elem):
        
            #get connectivity of current element
            elem_connect_i = tria_elem_connect[elem_ind,:]

            #get coordinate for current element
            coord_data_x = nc.variables['coordx'][elem_connect_i]
            coord_data_y = nc.variables['coordy'][elem_connect_i]

            #count (ptr of elem along czm)
            count = 0

            ##first ptr
            x1 = coord_data_x[0]; y1 = coord_data_y[0]

            #check whether this ptr belongs to czm
            find_pos_ind_y1 = np.where(abs(y1 - 0.0)<1)[0]
            if find_pos_ind_y1.size > 0:
                count += 1
        
            ##second ptr
            x2 = coord_data_x[1]; y2 = coord_data_y[1]

            #check whether this ptr belongs to czm
            find_pos_ind_y2 = np.where(abs(y2 - 0.0)<1)[0]
            if find_pos_ind_y2.size > 0:
                count += 1

            ##third ptr
            x3 = coord_data_x[2]; y3 = coord_data_y[2]

            #check whether this ptr belongs to czm
            find_pos_ind_y3 = np.where(abs(y3 - 0.0)<1)[0]
            if find_pos_ind_y3.size > 0:
                count += 1
        
            if count > 1: #this is an element along czm
                list_mf_elem.append(elem_ind)
                list_mf_xcoord.extend([x1,x2,x3])
                list_mf_ycoord.extend([y1,y2,y3])
                list_mf_xcoordcenter.extend([(x1+x2+x3)/3])
                list_mf_ycoordcenter.extend([(y1+y2+y3)/3])
        
        #save list_czm_elem
        np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(elemblock_ind)+".txt",list_mf_elem,fmt='%i',newline=" ")
        np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_xcoord"+str(elemblock_ind)+".txt",list_mf_xcoord,fmt='%i',newline=" ")
        np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_ycoord"+str(elemblock_ind)+".txt",list_mf_ycoord,fmt='%i',newline=" ")
        np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_xcoordcenter"+str(elemblock_ind)+".txt",list_mf_xcoordcenter,fmt='%i',newline=" ")
        np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_ycoordcenter"+str(elemblock_ind)+".txt",list_mf_ycoordcenter,fmt='%i',newline=" ")

    #check
    #subdomainIDcheck = len(list_czm_elem) * [100]
    #np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/list_czm_ids_check.txt", subdomainIDcheck, fmt='%i',newline=" ")

    #check
    plt.figure()
    for block in range(1,num_elem_block+1):
        if os.path.getsize("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(block)+".txt") != 0:
            list_mf_elem = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(block)+".txt")
            list_mf_xcoord = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_xcoord"+str(block)+".txt")
            list_mf_ycoord = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_ycoord"+str(block)+".txt")
            for elem in range(len(list_mf_elem)):
                plt.plot([list_mf_xcoord[3*elem],list_mf_xcoord[3*elem+1]],[list_mf_ycoord[3*elem],list_mf_ycoord[3*elem+1]])
                plt.plot([list_mf_xcoord[3*elem+1],list_mf_xcoord[3*elem+2]],[list_mf_ycoord[3*elem+1],list_mf_ycoord[3*elem+2]])
                plt.plot([list_mf_xcoord[3*elem+2],list_mf_xcoord[3*elem]],[list_mf_ycoord[3*elem+2],list_mf_ycoord[3*elem]])

    plt.show()


def RetrieveCZMElementIDs(exodus_file_path,
                          save_list_czm_elem_file_path,
                          dict_ptrs_data_file_path,
                          master_locs,
                          fault_name):

    """
    Retrieve all Element IDs for each CZM region

    note: 
    master_locs: prior knowledge of master surface location [num1 num2 num3 num4]
    such that " num1 * x_coord + num2 * y_coord + num3 (num4 0:> 1:<) 0 "
    """

    ###--------------------------Read czm_ptrs_data.txt--------------------------###
    #load txt
    data_ptrs = np.loadtxt(dict_ptrs_data_file_path, comments='"Points:0"', delimiter=',')

    #get column functions
    #['"Points:0"' '"Points:1"' '"Points:2"']
    data_ptrs_func = np.loadtxt(dict_ptrs_data_file_path, dtype='str', max_rows=1, delimiter=',')

    #get which col "Points:0"
    col_Points0 = np.where(data_ptrs_func == '"Points:0"')

    #get which col "Points:1"
    col_Points1 = np.where(data_ptrs_func == '"Points:1"')

    #get ptrs along on the branch faults
    data_ptrsxy = np.hstack((data_ptrs[:,col_Points0],data_ptrs[:,col_Points1]))

    ## check
    # plt.figure()
    # plt.plot(data_ptrsxy[:,0],data_ptrsxy[:,1],'.')
    # plt.show()
    # exit(0)

    ###--------------------------------------------------------------------------###

    ###--------------------------Read exodus file-----------------------------------###
    #read file
    nc = netCDF4.Dataset(exodus_file_path)

    #get number of element block
    num_elem_block = len(nc.variables['eb_status'])
    
    ###--------------------------------------------------------------------------###

    #loop over element block
    for elemblock_ind in range(1,num_elem_block+1): #don't skip the background elem ind 0

        #initilize list for saving elem along czm
        list_czm_elem = []

        #initialize list for saving x coordinate
        list_czm_elem_xcoord = []

        #initialize list for saving y coordinate
        list_czm_elem_ycoord = []

        #initialize list for saving centroid point x coord of tria element (for ploting)
        list_czm_elem_xcoordcenter = []

        #initialize list for saving centroid point y coord of tria element (for ploting)
        list_czm_elem_ycoordcenter = []
        
        #track progress
        print("(Main) elemblock_ind: ", elemblock_ind)
        print("(Main) Progress: ", elemblock_ind/(num_elem_block)*100,"%")

        #get connectivity for the current block
        #shrift by 1
        tria_elem_connect = nc.variables['connect' + str(elemblock_ind)][:]-1

        #TRIA3 Elem Num
        num_elem = np.shape(tria_elem_connect)[0]

        #loop over element
        for elem_ind in range(num_elem):
        
            #get connectivity of current element
            elem_connect_i = tria_elem_connect[elem_ind,:]

            #get coordinate for current element
            coord_data_x = nc.variables['coordx'][elem_connect_i]
            coord_data_y = nc.variables['coordy'][elem_connect_i]

            #count (ptr of elem along czm)
            count = 0

            ##first ptr
            x1 = coord_data_x[0]; y1 = coord_data_y[0]

            ##tol
            tol = 5.0 #1.5

            #check whether this ptr belongs to czm
            find_pos_ind_x1 = np.where(abs(x1 - data_ptrsxy[:,0])<tol)[0]
            find_pos_ind_y1 = np.where(abs(y1 - data_ptrsxy[:,1])<tol)[0]
            common_flag_ptr_1 = np.intersect1d(find_pos_ind_x1,find_pos_ind_y1)
            if common_flag_ptr_1.size > 0:
                count += 1
        
            ##second ptr
            x2 = coord_data_x[1]; y2 = coord_data_y[1]

            #check whether this ptr belongs to czm
            find_pos_ind_x2 = np.where(abs(x2 - data_ptrsxy[:,0])<tol)[0]
            find_pos_ind_y2 = np.where(abs(y2 - data_ptrsxy[:,1])<tol)[0]
            common_flag_ptr_2 = np.intersect1d(find_pos_ind_x2,find_pos_ind_y2)
            if common_flag_ptr_2.size > 0:
                count += 1

            ##third ptr
            x3 = coord_data_x[2]; y3 = coord_data_y[2]

            #check whether this ptr belongs to czm
            find_pos_ind_x3 = np.where(abs(x3 - data_ptrsxy[:,0])<tol)[0]
            find_pos_ind_y3 = np.where(abs(y3 - data_ptrsxy[:,1])<tol)[0]
            common_flag_ptr_3 = np.intersect1d(find_pos_ind_x3,find_pos_ind_y3)
            if common_flag_ptr_3.size > 0:
                count += 1

            if count > 1: #this is an element along czm
                
                #master_locs constraint
                elemcenterx = (x1+x2+x3)/3
                elemcentery = (y1+y2+y3)/3
                if master_locs[3] == 0: #>
                    elembelongstomaster = master_locs[0] * elemcenterx + master_locs[1] * elemcentery + master_locs[2] > 0
                elif master_locs[3] == 1: #<
                    elembelongstomaster = master_locs[0] * elemcenterx + master_locs[1] * elemcentery + master_locs[2] < 0

                if elembelongstomaster == True: #this is a master element along czm
            
                    list_czm_elem.append(elem_ind)
                    list_czm_elem_xcoord.extend([x1,x2,x3])
                    list_czm_elem_ycoord.extend([y1,y2,y3])
                    list_czm_elem_xcoordcenter.extend([(x1+x2+x3)/3])
                    list_czm_elem_ycoordcenter.extend([(y1+y2+y3)/3])
    
        #save list_czm_elem
        np.savetxt(save_list_czm_elem_file_path + "/list_czm_elem"+str(elemblock_ind)+".txt",list_czm_elem,fmt='%i',newline=" ")
        np.savetxt(save_list_czm_elem_file_path + "/list_czm_elem_xcoord"+str(elemblock_ind)+".txt",list_czm_elem_xcoord,fmt='%i',newline=" ")
        np.savetxt(save_list_czm_elem_file_path + "/list_czm_elem_ycoord"+str(elemblock_ind)+".txt",list_czm_elem_ycoord,fmt='%i',newline=" ")
        np.savetxt(save_list_czm_elem_file_path + "/list_czm_elem_xcoordcenter"+str(elemblock_ind)+".txt",list_czm_elem_xcoordcenter,fmt='%i',newline=" ")
        np.savetxt(save_list_czm_elem_file_path + "/list_czm_elem_ycoordcenter"+str(elemblock_ind)+".txt",list_czm_elem_ycoordcenter,fmt='%i',newline=" ")

    #check
    #subdomainIDcheck = len(list_czm_elem) * [100]
    #np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/list_czm_ids_check.txt", subdomainIDcheck, fmt='%i',newline=" ")

    #check
    plt.figure()
    for block in range(1,num_elem_block+1): #don't remove background elem ind 0
        if os.path.getsize(save_list_czm_elem_file_path + "/list_czm_elem"+str(block)+".txt") != 0:
            list_czm_elem = np.loadtxt(save_list_czm_elem_file_path + "/list_czm_elem"+str(block)+".txt")
            list_czm_elem_xcoord = np.loadtxt(save_list_czm_elem_file_path + "/list_czm_elem_xcoord"+str(block)+".txt")
            list_czm_elem_ycoord = np.loadtxt(save_list_czm_elem_file_path + "/list_czm_elem_ycoord"+str(block)+".txt")
            #hardcode #the problem here is branch point is also here
            if fault_name == "mf180" and block == 46:
                continue
            elif fault_name == "mf30" and ( block == 40 or block == 41 or block == 42 or block == 44 ):
                continue
            else:
                # print(block)
                for elem in range(len(list_czm_elem)):
                    plt.plot([list_czm_elem_xcoord[3*elem],list_czm_elem_xcoord[3*elem+1]],[list_czm_elem_ycoord[3*elem],list_czm_elem_ycoord[3*elem+1]])
                    plt.plot([list_czm_elem_xcoord[3*elem+1],list_czm_elem_xcoord[3*elem+2]],[list_czm_elem_ycoord[3*elem+1],list_czm_elem_ycoord[3*elem+2]])
                    plt.plot([list_czm_elem_xcoord[3*elem+2],list_czm_elem_xcoord[3*elem]],[list_czm_elem_ycoord[3*elem+2],list_czm_elem_ycoord[3*elem]])

    plt.show()


def LinePlotParamVal3D(exodus_file_path,plot_var_name,plot_var_threshold,plot_flag):

    """
    3D Line Plot Slip Rate Plot
    """

    ##
    #Initialize dict storing elem var name / index pair
    dict_evn_index = {}

    #Read Elem Var Name
    evn = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_name_elem_var.txt",dtype=str)

    #fill in dict_evn_index
    for evn_ind in range(len(evn)):
        evn_i = evn[evn_ind]
        # start from 1
        dict_evn_index[evn_i] = evn_ind + 1

    #get required param to be plotted
    var_ind = dict_evn_index[plot_var_name]

    ##
    #Read Exodus File
    nc = netCDF4.Dataset(exodus_file_path)

    #Obtain Time Series
    timeseries = nc.variables['time_whole']

    ##
    #get number of element block
    num_elem_block = len(nc.variables['eb_status'])

    #czm#
    #initialize dict storing param time series array
    dict_blockid_arrparamval = {}

    #initialize dict storing elem xcoordcenter
    dict_blockid_xcoordcenter = {}

    #initialize dict storing elem ycoordcenter
    dict_blockid_ycoordcenter = {}

    #mf#
    #initialize dict storing param time series array
    dict_blockid_arrparamval_mf = {}

    #initialize dict storing elem xcoordcenter
    dict_blockid_xcoordcenter_mf = {}

    #initialize dict storing elem ycoordcenter
    dict_blockid_ycoordcenter_mf = {}

    #czm
    if plot_flag == "czmonly":
        
        #Loop Over Time
        for time_ind in range(len(timeseries)):
        #for time_ind in [0, 200]:

            #Get current time
            time_i = timeseries[time_ind]

            #round 2 decimal
            time_i = np.round(time_i,2)

            #Create figure
            # fig = plt.figure(figsize=(100,10),constrained_layout=True)
            # ax = fig.add_subplot(projection = '3d')
            plt.figure(figsize=(50,10))

            #initialize arr storing all data
            list_czm_elem_xcoordcenter_i_all = []
            list_czm_elem_ycoordcenter_i_all = []
            list_czm_param_i_all = []

            #Loop Over Blocks
            for elemblock_ind in range(1,num_elem_block+1):

                #track
                print("[Time]: ", time_i)
                print("(Sub) elemblock_ind: ", elemblock_ind)
                print("(Sub) Progress: ", elemblock_ind/(num_elem_block)*100,"%")

                #load/store in dict at initial time
                if time_i == 0.0:

                    #load elem ind list
                    list_czm_elem_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_czm_elem"+str(elemblock_ind)+".txt")

                    #int
                    list_czm_elem_i = [int(czm_elem_i) for czm_elem_i in list_czm_elem_i]

                    #load elem center point coordinate
                    list_czm_elem_xcoordcenter_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_czm_elem_xcoordcenter"+str(elemblock_ind)+".txt")
                    list_czm_elem_ycoordcenter_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_czm_elem_ycoordcenter"+str(elemblock_ind)+".txt")

                    #Obtain required param store array
                    arr_param_i = nc.variables['vals_elem_var'+str(var_ind)+'eb'+str(elemblock_ind)]

                    #Extract needed values from arr_param_i (with time series) (time, element_i value)
                    arr_czm_param_all = arr_param_i[:,list_czm_elem_i]

                    dict_blockid_xcoordcenter[elemblock_ind] = np.asarray(list_czm_elem_xcoordcenter_i)
                    dict_blockid_ycoordcenter[elemblock_ind] = np.asarray(list_czm_elem_ycoordcenter_i)
                    dict_blockid_arrparamval[elemblock_ind]  = arr_czm_param_all

                    #Get param val for current time step
                    arr_czm_param_i = arr_czm_param_all[time_ind, :]

                #load existing variables in dict
                else:

                    #load xcoordcenter
                    list_czm_elem_xcoordcenter_i = dict_blockid_xcoordcenter[elemblock_ind]

                    #load ycoordcenter
                    list_czm_elem_ycoordcenter_i = dict_blockid_ycoordcenter[elemblock_ind]

                    #load arrparamval
                    arr_czm_param_i = dict_blockid_arrparamval[elemblock_ind][time_ind, :]

                #find param less than threshold, treat as zero
                ind_threshold_pos = np.argwhere( (arr_czm_param_i >  plot_var_threshold) )
                ind_threshold_neg = np.argwhere( (arr_czm_param_i < -plot_var_threshold) )

                ind_threshold = np.vstack((ind_threshold_pos,ind_threshold_neg))

                #update list/arr
                list_czm_elem_xcoordcenter_i = list_czm_elem_xcoordcenter_i[ind_threshold]
                list_czm_elem_ycoordcenter_i = list_czm_elem_ycoordcenter_i[ind_threshold]
                list_czm_param_i = arr_czm_param_i[ind_threshold].tolist()

                #store in temp list
                list_czm_elem_xcoordcenter_i_all.extend(list_czm_elem_xcoordcenter_i)
                list_czm_elem_ycoordcenter_i_all.extend(list_czm_elem_ycoordcenter_i)
                list_czm_param_i_all.extend(list_czm_param_i)

            #plot line
            my_cmap = plt.get_cmap('jet')
            p = plt.scatter(list_czm_elem_xcoordcenter_i_all,list_czm_elem_ycoordcenter_i_all,s=20.0, c=list_czm_param_i_all, cmap = my_cmap, vmin= -0.2, vmax=0.6)
            plt.colorbar(p)
            plt.grid(False)
            plt.xlim([-10000,10000])
            plt.ylim([-1000,1000])
            plt.xlabel("x coord")
            plt.ylabel("y coord")
            plt.title("Branch Faults Slip Rate Spatial Distribution at Time Step "+str(time_i)+" s") #name needed to be generalized
            plt.savefig("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output_png/sliprate_czmonly_"+str(time_ind)+".png",bbox_inches = 'tight')
            #plt.show()
            #exit(0)
    
    #both branch faults and main faults
    elif plot_flag == "mfczm":

        #Loop Over Time
        for time_ind in range(len(timeseries)):
        #for time_ind in [0, 200]:

            #Get current time
            time_i = timeseries[time_ind]

            #round 2 decimal
            time_i = np.round(time_i,2)

            #Create figure
            # fig = plt.figure(figsize=(100,10),constrained_layout=True)
            # ax = fig.add_subplot(projection = '3d')
            plt.figure(figsize=(50,10))

            #initialize arr storing all data
            #czm + mf#
            list_elem_xcoordcenter_i_all = []
            list_elem_ycoordcenter_i_all = []
            list_param_i_all = []

            #Loop Over Blocks
            for elemblock_ind in range(1,num_elem_block+1):

                #track
                print("[Time]: ", time_i)
                print("(Sub) elemblock_ind: ", elemblock_ind)
                print("(Sub) Progress: ", elemblock_ind/(num_elem_block)*100,"%")

                #load/store in dict at initial time
                if time_i == 0.0:
                    
                    #czm#
                    #load elem ind list
                    list_czm_elem_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_czm_elem"+str(elemblock_ind)+".txt")

                    #int
                    list_czm_elem_i = [int(czm_elem_i) for czm_elem_i in list_czm_elem_i]

                    #load elem center point coordinate
                    list_czm_elem_xcoordcenter_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_czm_elem_xcoordcenter"+str(elemblock_ind)+".txt")
                    list_czm_elem_ycoordcenter_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_czm_elem_ycoordcenter"+str(elemblock_ind)+".txt")
                    ##

                    #mf#
                    if os.path.getsize("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(elemblock_ind)+".txt") != 0:
                        #load elem ind list
                        list_mf_elem_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(elemblock_ind)+".txt")

                        #int
                        list_mf_elem_i = [int(mf_elem_i) for mf_elem_i in list_mf_elem_i]

                        #load elem center point coordinate
                        list_mf_elem_xcoordcenter_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_xcoordcenter"+str(elemblock_ind)+".txt")
                        list_mf_elem_ycoordcenter_i = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_ycoordcenter"+str(elemblock_ind)+".txt")
                    ##    

                    #Obtain required param store array
                    arr_param_i = nc.variables['vals_elem_var'+str(var_ind)+'eb'+str(elemblock_ind)]
                    
                    #czm#
                    #Extract needed values from arr_param_i (with time series) (time, element_i value)
                    arr_czm_param_all = arr_param_i[:,list_czm_elem_i]

                    dict_blockid_xcoordcenter[elemblock_ind] = np.asarray(list_czm_elem_xcoordcenter_i)
                    dict_blockid_ycoordcenter[elemblock_ind] = np.asarray(list_czm_elem_ycoordcenter_i)
                    dict_blockid_arrparamval[elemblock_ind]  = arr_czm_param_all

                    #Get param val for current time step
                    arr_czm_param_i = arr_czm_param_all[time_ind, :]

                    #mf#
                    if os.path.getsize("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(elemblock_ind)+".txt") != 0:
                        #Extract needed values from arr_param_i (with time series) (time, element_i value)
                        arr_mf_param_all = arr_param_i[:,list_mf_elem_i]

                        #attach to existing array in dict
                        dict_blockid_xcoordcenter_mf[elemblock_ind] = np.asarray(list_mf_elem_xcoordcenter_i)
                        dict_blockid_ycoordcenter_mf[elemblock_ind] = np.asarray(list_mf_elem_ycoordcenter_i)
                        dict_blockid_arrparamval_mf[elemblock_ind]  = arr_mf_param_all

                        #Get param val for current time step
                        arr_mf_param_i = arr_mf_param_all[time_ind, :]

                #load existing variables in dict
                else:
                    
                    #czm#
                    #load xcoordcenter
                    list_czm_elem_xcoordcenter_i = dict_blockid_xcoordcenter[elemblock_ind]

                    #load ycoordcenter
                    list_czm_elem_ycoordcenter_i = dict_blockid_ycoordcenter[elemblock_ind]

                    #load arrparamval
                    arr_czm_param_i = dict_blockid_arrparamval[elemblock_ind][time_ind, :]

                    #mf#
                    if os.path.getsize("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(elemblock_ind)+".txt") != 0:
                        #load xcoordcenter
                        list_mf_elem_xcoordcenter_i = dict_blockid_xcoordcenter_mf[elemblock_ind]

                        #load ycoordcenter
                        list_mf_elem_ycoordcenter_i = dict_blockid_ycoordcenter_mf[elemblock_ind]

                        #load arrparamval
                        arr_mf_param_i = dict_blockid_arrparamval_mf[elemblock_ind][time_ind, :]

                #czm#
                #find param less than threshold, treat as zero
                ind_threshold_pos = np.argwhere( (arr_czm_param_i >  plot_var_threshold) )
                ind_threshold_neg = np.argwhere( (arr_czm_param_i < -plot_var_threshold) )

                ind_threshold = np.vstack((ind_threshold_pos,ind_threshold_neg))

                #update list/arr
                list_czm_elem_xcoordcenter_i = list_czm_elem_xcoordcenter_i[ind_threshold]
                list_czm_elem_ycoordcenter_i = list_czm_elem_ycoordcenter_i[ind_threshold]
                list_czm_param_i = arr_czm_param_i[ind_threshold].tolist()

                #store in temp list
                list_elem_xcoordcenter_i_all.extend(list_czm_elem_xcoordcenter_i)
                list_elem_ycoordcenter_i_all.extend(list_czm_elem_ycoordcenter_i)
                list_param_i_all.extend(list_czm_param_i)

                #mf#
                if os.path.getsize("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output/list_mf_elem"+str(elemblock_ind)+".txt") != 0:
                    #find param less than threshold, treat as zero
                    ind_threshold_pos = np.argwhere( (arr_mf_param_i >  plot_var_threshold) )
                    ind_threshold_neg = np.argwhere( (arr_mf_param_i < -plot_var_threshold) )

                    ind_threshold = np.vstack((ind_threshold_pos,ind_threshold_neg))

                    #update list/arr
                    list_mf_elem_xcoordcenter_i = list_mf_elem_xcoordcenter_i[ind_threshold]
                    list_mf_elem_ycoordcenter_i = list_mf_elem_ycoordcenter_i[ind_threshold]
                    list_mf_param_i = arr_mf_param_i[ind_threshold].tolist()

                    #store in temp list
                    list_elem_xcoordcenter_i_all.extend(list_mf_elem_xcoordcenter_i)
                    list_elem_ycoordcenter_i_all.extend(list_mf_elem_ycoordcenter_i)
                    list_param_i_all.extend(list_mf_param_i)

            #plot line
            my_cmap = plt.get_cmap('jet')
            p = plt.scatter(list_elem_xcoordcenter_i_all,list_elem_ycoordcenter_i_all,s=20.0, c=list_param_i_all, cmap = my_cmap, vmin= -1.0, vmax=4.0)
            plt.colorbar(p)
            plt.grid(False)
            plt.xlim([-10000,10000])
            plt.ylim([-1000,1000])
            plt.xlabel("x coord")
            plt.ylabel("y coord")
            plt.title("Slip Rate Spatial Distribution at Time Step "+str(time_i)+" s") #name needed to be generalized
            plt.savefig("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv5_contmfv2/postv5cmfv2/folder_output_png/sliprate_mfczm_"+str(time_ind)+".png",bbox_inches = 'tight')
            #plt.show()
            #exit(0)

def LinePlotParamVal2D(exodus_file_path,
                       plot_var_name,
                       save_folder_output_file_path,
                       save_folder_output_png_file_path,
                       fault_name,
                       additional_files = None,
                       additional_flags = None,
                       dict_ini_shear_sts = None,
                       dict_ini_normal_sts = None):

    """
    2D Line Plot Variables

    Note:
    fault_name:
        *current on which fault [mf180, mf30, bf30]
    additional_files: 
        *front_position [rupture front position vs time](need to choose plot_var_name = 'tangent_jump_rate')
            {mf180/mf30/bf30: [(rupture starts positon x, rupture starts positon y), (rupture end position x, rupture end position y)]}
    additional_flags:
        *front_position (string name)
    dict_ini_shear_sts: initial shear stress for each fault
    """

    ##
    #Initialize dict storing elem var name / index pair
    dict_evn_index = {}

    #Read Elem Var Name
    evn = np.loadtxt(save_folder_output_file_path + "/list_name_elem_var.txt",dtype=str)

    #fill in dict_evn_index
    for evn_ind in range(len(evn)):
        evn_i = evn[evn_ind]
        # start from 1
        dict_evn_index[evn_i] = evn_ind + 1

    #current fault name
    print("[current fault name] :", fault_name)

    #get required param to be plotted
    #loop over plot names
    for plot_var_name_index in range(len(plot_var_name)):
        
        plot_var_name_i = plot_var_name[plot_var_name_index]
        var_ind = dict_evn_index[plot_var_name_i]

        ##
        #Read Exodus File
        nc = netCDF4.Dataset(exodus_file_path)

        #Obtain Time Series
        timeseries = nc.variables['time_whole']

        ##
        #get number of element block
        num_elem_block = len(nc.variables['eb_status'])

        #czm#
        #initialize dict storing param time series array
        dict_blockid_arrparamval = {}

        #other dir if needed
        dict_blockid_arrparamval_other = {}

        #initialize dict storing elem xcoordcenter
        dict_blockid_xcoordcenter = {}

        #initialize dict storing elem ycoordcenter
        dict_blockid_ycoordcenter = {}

        #For front_position/tangent_jump_rate
        if (plot_var_name_i == "tangent_jump_rate" and additional_flags == "front_position") or (plot_var_name_i == "jump_x_rate" and additional_flags == "front_position"):

            #get number of nodes
            numofnodes = len(additional_files[fault_name])

            print(numofnodes)
            
            #initialize node_pair_list
            node_pair_list_x = []
            node_pair_list_y = []

            #initialize front-time list
            if numofnodes == 4:
                list_front_time0 = []
                list_front_xcoord0 = []
                list_front_time1 = []
                list_front_xcoord1 = []
            else:
                list_front_time0 = []
                list_front_xcoord0 = []

            #get node coordinate pairs
            count = 0
            for node_pair_i in range(int(numofnodes/2)):

                node_k = additional_files[fault_name][node_pair_i + count]
                node_l = additional_files[fault_name][node_pair_i + count + 1]

                node_pair_list_x.append([node_k[0],node_l[0]])
                node_pair_list_y.append([node_k[1],node_l[1]])

                count += 1
            

            # print(save_folder_output_file_path)

        #For global jump_rate, jump
        if plot_var_name_i == "jump_x_rate":
            jump_y_rate_ind = dict_evn_index["jump_y_rate"]
        if plot_var_name_i == "jump_x":
            jump_y_ind = dict_evn_index["jump_y"]
        if plot_var_name_i == "traction_x":
            traction_y_ind = dict_evn_index["traction_y"]

        print("[current plot name]", plot_var_name_i)

        #record maximum slip and average slip
        if plot_var_name_i == "tangent_jump":

            list_slip = []

        #Loop Over Time
        for time_ind in range(len(timeseries)):
        # for time_ind in [0]:

            #Get current time
            time_i = timeseries[time_ind]

            #round 2 decimal
            time_i = np.round(time_i,2)

            print("[Time]: ", time_i)

            #Create figure
            # fig = plt.figure(figsize=(100,10),constrained_layout=True)
            # ax = fig.add_subplot(projection = '3d')
            plt.figure(figsize=(30,10))

            #initialize arr storing all data
            list_czm_elem_xcoordcenter_i_all = []
            list_czm_elem_ycoordcenter_i_all = []
            list_czm_param_i_all = []

            #other
            list_czm_param_j_all = []

            #Loop Over Blocks
            #don't remove background elem ind 0
            for elemblock_ind in range(1,num_elem_block+1):

                #track
                # print("[Time]: ", time_i)
                # print("(Sub) elemblock_ind: ", elemblock_ind)
                # print("(Sub) Progress: ", elemblock_ind/(num_elem_block)*100,"%")

                #check if it is empty
                #load elem ind list #hardcode
                list_czm_elem_i = np.loadtxt(save_folder_output_file_path + "/list_czm_elem"+str(elemblock_ind)+".txt")
                if elemblock_ind == 46 and fault_name == "mf180":
                    continue
                elif fault_name == "mf30" and ( elemblock_ind == 40 or elemblock_ind == 41 or elemblock_ind == 42 or elemblock_ind == 44 ):
                    continue
                elif len(list_czm_elem_i) == 0: #empty list, abadon this block
                    continue
                else:
                    #load/store in dict at initial time
                    if time_i == 0.0:

                        #load elem ind list
                        list_czm_elem_i = np.loadtxt(save_folder_output_file_path + "/list_czm_elem"+str(elemblock_ind)+".txt")

                        #int
                        list_czm_elem_i = [int(czm_elem_i) for czm_elem_i in list_czm_elem_i]

                        #load elem center point coordinate
                        list_czm_elem_xcoordcenter_i = np.loadtxt(save_folder_output_file_path + "/list_czm_elem_xcoordcenter"+str(elemblock_ind)+".txt")
                        print(elemblock_ind)
                        list_czm_elem_ycoordcenter_i = np.loadtxt(save_folder_output_file_path + "/list_czm_elem_ycoordcenter"+str(elemblock_ind)+".txt")

                        #Obtain required param store array
                        arr_param_i = nc.variables['vals_elem_var'+str(var_ind)+'eb'+str(elemblock_ind)]

                        #global-local jump, jump_rate, traction
                        if fault_name != "mf180":
                            if plot_var_name_i == "jump_x_rate" or plot_var_name_i == "traction_x" or plot_var_name_i == "jump_x":

                                arr_param_i = arr_param_i[:]

                                if plot_var_name_i == "jump_x_rate":
                                    arr_param_j = nc.variables['vals_elem_var'+str(jump_y_rate_ind)+'eb'+str(elemblock_ind)]
                                if plot_var_name_i == "traction_x":
                                    arr_param_j = nc.variables['vals_elem_var'+str(traction_y_ind)+'eb'+str(elemblock_ind)]
                                if plot_var_name_i == "jump_x":
                                    arr_param_j = nc.variables['vals_elem_var'+str(jump_y_ind)+'eb'+str(elemblock_ind)]
                                
                                for ind in range(len(arr_param_i)):

                                    gparam_x_i = arr_param_i[ind]
                                    gparam_y_i = arr_param_j[ind]

                                    if fault_name == "bf30":

                                        cost = math.cos(math.radians(32))
                                        sint = math.sin(math.radians(32))
                                    
                                        arr_param_i[ind] = gparam_x_i * cost + gparam_y_i * sint

                                    elif fault_name == "mf30":

                                        cost = math.cos(math.radians(35))
                                        sint = math.sin(math.radians(35))

                                        arr_param_i[ind] =  gparam_x_i * cost + gparam_y_i * sint
                        else:
                            if plot_var_name_i == "traction_x":
                                    arr_param_j = nc.variables['vals_elem_var'+str(traction_y_ind)+'eb'+str(elemblock_ind)]
                                    
                                    
                        #Extract needed values from arr_param_i (with time series) (time, element_i value)
                        arr_czm_param_all = arr_param_i[:,list_czm_elem_i]

                        dict_blockid_xcoordcenter[elemblock_ind] = np.asarray(list_czm_elem_xcoordcenter_i)
                        dict_blockid_ycoordcenter[elemblock_ind] = np.asarray(list_czm_elem_ycoordcenter_i)
                        dict_blockid_arrparamval[elemblock_ind]  = arr_czm_param_all

                        #Get param val for current time step
                        arr_czm_param_i = arr_czm_param_all[time_ind, :]

                        #y-dir if needed
                        if plot_var_name_i == "traction_x" and additional_flags == "stress_and_strength":
                            arr_czm_param_all_other = arr_param_j[:,list_czm_elem_i]
                            dict_blockid_arrparamval_other[elemblock_ind] = arr_czm_param_all_other
                            arr_czm_param_j = arr_czm_param_all[time_ind, :]


                    #load existing variables in dict
                    else:

                        #load xcoordcenter
                        list_czm_elem_xcoordcenter_i = dict_blockid_xcoordcenter[elemblock_ind]

                        #load ycoordcenter
                        list_czm_elem_ycoordcenter_i = dict_blockid_ycoordcenter[elemblock_ind]

                        #load arrparamval
                        arr_czm_param_i = dict_blockid_arrparamval[elemblock_ind][time_ind, :]

                        #load other dir
                        if plot_var_name_i == "traction_x" and additional_flags == "stress_and_strength":
                            arr_czm_param_j = dict_blockid_arrparamval_other[elemblock_ind][time_ind, :]

                    #update list/arr
                    list_czm_param_i = arr_czm_param_i.tolist()

                    #store in temp list
                    list_czm_elem_xcoordcenter_i_all.extend(list_czm_elem_xcoordcenter_i)
                    list_czm_elem_ycoordcenter_i_all.extend(list_czm_elem_ycoordcenter_i)
                    list_czm_param_i_all.extend(list_czm_param_i)

                    #update/store other
                    if plot_var_name_i == "traction_x" and additional_flags == "stress_and_strength":
                        list_czm_param_j = arr_czm_param_j.tolist()
                        list_czm_param_j_all.extend(list_czm_param_j)

                #sort array
                sort_index = np.argsort(list_czm_elem_xcoordcenter_i_all)
                sorted_arr_czm_elem_xcoordcenter_i_all = np.array(list_czm_elem_xcoordcenter_i_all)[sort_index]
                sorted_arr_czm_elem_ycoordcenter_i_all = np.array(list_czm_elem_ycoordcenter_i_all)[sort_index]

                #compute arc-length
                if fault_name == "bf30":
                    sorted_arr_czm_elem_arclength_i_all = np.zeros(np.shape(sorted_arr_czm_elem_xcoordcenter_i_all))
                    for ind in range(np.shape(sorted_arr_czm_elem_xcoordcenter_i_all)[0]):
                        xcoord_i = sorted_arr_czm_elem_xcoordcenter_i_all[ind]
                        ycoord_i = sorted_arr_czm_elem_ycoordcenter_i_all[ind]
                        arclen_i = np.sqrt(xcoord_i ** 2 + ycoord_i ** 2)
                        if xcoord_i < 0:
                            arclen_i = - arclen_i
                        sorted_arr_czm_elem_arclength_i_all[ind] = arclen_i
                elif fault_name == "mf180":
                    sorted_arr_czm_elem_arclength_i_all = sorted_arr_czm_elem_xcoordcenter_i_all
                elif fault_name == "mf30":
                    endptr_x = additional_files[fault_name][0][0]
                    endptr_y = additional_files[fault_name][0][1]
                    sorted_arr_czm_elem_arclength_i_all = np.zeros(np.shape(sorted_arr_czm_elem_xcoordcenter_i_all))
                    for ind in range(np.shape(sorted_arr_czm_elem_xcoordcenter_i_all)[0]):
                        xcoord_i = sorted_arr_czm_elem_xcoordcenter_i_all[ind]
                        ycoord_i = sorted_arr_czm_elem_ycoordcenter_i_all[ind]
                        arclen_i = np.sqrt( (endptr_x - xcoord_i) ** 2 + (endptr_y - ycoord_i) ** 2)
                        if xcoord_i < 0:
                            arclen_i = - arclen_i
                        sorted_arr_czm_elem_arclength_i_all[ind] = arclen_i

                sorted_arr_czm_param_i_all = np.array(list_czm_param_i_all)[sort_index]
                # print(sorted_arr_czm_param_i_all)

                #sort other
                if plot_var_name_i == "traction_x" and additional_flags == "stress_and_strength":
                    sorted_arr_czm_param_j_all = np.array(list_czm_param_j_all)[sort_index]

                if plot_var_name_i == "tangent_jump":

                    #output maximum and average
                    print("maximum slip: ", np.max(sorted_arr_czm_param_i_all))
                    print("average slip: ", np.average(sorted_arr_czm_param_i_all))

                    list_slip.append(np.max(sorted_arr_czm_param_i_all)) 
                    list_slip.append(np.average(sorted_arr_czm_param_i_all))

                #plot line
                #tangent dir
                if plot_var_name_i == "tangent_jump_rate" and additional_flags == None:
                    # plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all,'b*-',linewidth=3.0)
                    # plt.grid(False)
                    # plt.xlabel("arc length (m)",fontsize=20)
                    # plt.ylim([-0.1,10.0])
                    # plt.ylabel("slip rate (m/s)",fontsize=20)
                    # if fault_name == "mf180":
                    #     plt.title("Main Fault Slip Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    # elif fault_name == "mf30":
                    #     plt.title("Fault Beyond Kink Slip Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    # else:
                    #     plt.title("Splay Fault Slip Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    # plt.savefig(save_folder_output_png_file_path + "/sliprate_"+str(time_ind)+".png",bbox_inches = 'tight')
                    # plt.xticks(fontsize=15)
                    # plt.yticks(fontsize=15)
                    # print("pass")
                    np.savetxt("./outputs/tangent_jump_rate/xcoord.txt",sorted_arr_czm_elem_arclength_i_all)
                    np.savetxt("./outputs/tangent_jump_rate/tangent_jump_rate_"+str(time_ind)+".txt",sorted_arr_czm_param_i_all)
                elif plot_var_name_i == "jump_x_rate" and additional_flags == None:
                    plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all,'b-',linewidth=3.0)
                    plt.grid(False)
                    plt.xlabel("arc length (m)",fontsize=20)
                    plt.ylim([-0.1,10.0])
                    plt.ylabel("slip rate (m/s)",fontsize=20)
                    if fault_name == "mf180":
                        plt.title("Main Fault Slip Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Slip Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Slip Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/sliprateglobaltolocal_"+str(time_ind)+".png",bbox_inches = 'tight')
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                elif plot_var_name_i == "tangent_jump":
                    # plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all,'b-',linewidth=3.0)
                    # plt.grid(False)
                    # plt.xlabel("arc length (m)",fontsize=20)
                    # plt.ylabel("slip (m)",fontsize=20)
                    # plt.ylim([-0.1,10.0])
                    # plt.xticks(fontsize=15)
                    # plt.yticks(fontsize=15)
                    # if fault_name == "mf180":
                    #     plt.title("Main Fault Slip at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    # elif fault_name == "mf30":
                    #     plt.title("Fault Beyond Kink Slip at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    # else:
                    #     plt.title("Splay Fault Slip at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    # # plt.savefig(save_folder_output_png_file_path + "/slip_"+str(time_ind)+".png")
                    np.savetxt("./outputs/tangent_jump/xcoord.txt",sorted_arr_czm_elem_arclength_i_all)
                    np.savetxt("./outputs/tangent_jump/tangent_jump_"+str(time_ind)+".txt",sorted_arr_czm_param_i_all)
                elif plot_var_name_i == "jump_x":
                    plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all,'b-',linewidth=3.0)
                    plt.grid(False)
                    plt.xlabel("arc length (m)",fontsize=20)
                    plt.ylabel("slip (m)",fontsize=20)
                    plt.ylim([-0.1,10.0])
                    if fault_name == "mf180":
                        plt.title("Main Fault Slip at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Slip at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Slip at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/slipg_"+str(time_ind)+".png",bbox_inches = 'tight')
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                elif plot_var_name_i == "tangent_traction":
                    if fault_name == "bf30":
                       #nucleation region 
                        nucleate_xstart = -26119.881 #hardcode
                        nucleate_ystart = -16321.513 #hardcode
                        nucleate_xend = -24763.005   #hardcode
                        nucleate_yend = -15473.643   #hardcode
                        nucleate_arclenstart = -1 * np.sqrt(nucleate_xstart**2 + nucleate_ystart**2) #hardcode
                        nucleate_arclenend = -1 * np.sqrt(nucleate_xend**2 + nucleate_yend**2)       #hardcode
                        nucleate_shear_sts = 24.7
                        #find index
                        range_index = np.where(np.logical_and(sorted_arr_czm_elem_arclength_i_all > nucleate_arclenstart, sorted_arr_czm_elem_arclength_i_all < nucleate_arclenend))[0]
                        #generate initial shear stress
                        list_ini_shear_sts = [dict_ini_shear_sts[fault_name]] * len(sorted_arr_czm_elem_arclength_i_all)
                        #fit in list with nucleate shear stress
                        for index in range_index:
                            list_ini_shear_sts[index] = nucleate_shear_sts
                        #plot
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,-sorted_arr_czm_param_i_all/1e6 + list_ini_shear_sts,'b-',linewidth=3.0)
                    else:
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,-sorted_arr_czm_param_i_all/1e6 + dict_ini_shear_sts[fault_name] ,'b-',linewidth=3.0)
                    plt.grid(False)
                    if fault_name == "mf180":
                        plt.xlim([-30000,50000])
                        plt.ylim([dict_ini_shear_sts[fault_name]-40,dict_ini_shear_sts[fault_name]+40])
                    elif fault_name == "mf30":
                        plt.xlim([-95532.164, -30000])
                        plt.ylim([dict_ini_shear_sts[fault_name]-40,dict_ini_shear_sts[fault_name]+40])
                    else:
                        plt.xlim([-33921.924,-160.033])
                        plt.ylim([dict_ini_shear_sts[fault_name]-20,dict_ini_shear_sts[fault_name]+20])
                    plt.xlabel("x coord (m)",fontsize=20)
                    plt.ylabel("Tragent Traction (MPa)",fontsize=20)
                    if fault_name == "mf180":
                        plt.title("Main Fault Tangent Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Tangent Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Tangent Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/tangtraction_"+str(time_ind)+".png",bbox_inches = 'tight')
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                elif plot_var_name_i == "traction_x" and additional_flags == None:
                    if fault_name == "bf30":
                        #nucleation region #hardcode #this is wrong, need care!!!
                        nucleate_xstart = -26119.881 #hardcode
                        nucleate_ystart = -16321.513 #hardcode
                        nucleate_xend = -24763.005   #hardcode
                        nucleate_yend = -15473.643   #hardcode
                        nucleate_arclenstart = -1 * np.sqrt(nucleate_xstart**2 + nucleate_ystart**2) #hardcode
                        nucleate_arclenend = -1 * np.sqrt(nucleate_xend**2 + nucleate_yend**2)       #hardcode
                        nucleate_shear_sts = 24.7
                        #find index
                        range_index = np.where(np.logical_and(sorted_arr_czm_elem_arclength_i_all > nucleate_arclenstart, sorted_arr_czm_elem_arclength_i_all < nucleate_arclenend))[0]
                        # print(range_index)
                        #generate initial shear stress
                        list_ini_shear_sts = [dict_ini_shear_sts[fault_name]] * len(sorted_arr_czm_elem_arclength_i_all)
                        # print(sorted_arr_czm_param_i_all[46:55]/1e6)
                        #fit in list with nucleate shear stress
                        # minindex = min(range_index); maxindex = max(range_index)
                        # print(minindex,maxindex)
                        # print(sorted_arr_czm_param_i_all[46:55]/1e6)
                        for index in [46,47,48,49,50,51,52,53,54]:
                            list_ini_shear_sts[index] = nucleate_shear_sts
                        #hardcode
                        list_ini_shear_sts[46] = nucleate_shear_sts - 3.198
                        list_ini_shear_sts[54] = nucleate_shear_sts - 3.198
                        # print(sorted_arr_czm_param_i_all[46:55]/1e6 + list_ini_shear_sts[46:55])
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all/1e6 + list_ini_shear_sts ,'b*-',linewidth=3.0)
                        # plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all/1e6 + dict_ini_shear_sts[fault_name] ,'b*-',linewidth=3.0)
                        plt.ylim([dict_ini_shear_sts[fault_name]-40,dict_ini_shear_sts[fault_name]+40])
                        plt.axhline(y=23.723)
                        plt.axhline(y=11.082)
                    elif fault_name == "mf180":
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all/1e6 + dict_ini_shear_sts[fault_name] ,'b*-',linewidth=3.0)
                        plt.ylim([dict_ini_shear_sts[fault_name]-40,dict_ini_shear_sts[fault_name]+40])
                    else:
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all/1e6 + dict_ini_shear_sts[fault_name] ,'b*-',linewidth=3.0)
                        plt.ylim([dict_ini_shear_sts[fault_name]-40,dict_ini_shear_sts[fault_name]+40])
                        plt.axhline(y=22.541)
                        plt.axhline(y=9.60)
                    plt.xlabel("Arc Length (m)",fontsize=20)
                    plt.ylabel("Tragent Traction (MPa)",fontsize=20)
                    if fault_name == "mf180":
                        plt.title("Main Fault Tangent Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Tangent Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Tangent Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/tangtraction_"+str(time_ind)+".png",bbox_inches = 'tight')
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                elif plot_var_name_i == "traction_x" and additional_flags == "stress_and_strength":
                    if fault_name == "mf180":
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all/1e6 + dict_ini_shear_sts[fault_name] ,'b*-',linewidth=3.0)
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,(sorted_arr_czm_param_j_all/1e6 + dict_ini_normal_sts[fault_name])*0.7 ,'r.-',linewidth=3.0)
                        plt.ylim([dict_ini_shear_sts[fault_name]-40,dict_ini_shear_sts[fault_name]+40])
                    if fault_name == "mf180":
                        plt.title("Main Fault Shear Stress / Shear Strength at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.legend(['Shear Stress', 'Shear Strength'],fontsize=20)
                    plt.xlabel("Arc Length (m)",fontsize=20)
                    plt.ylabel("Tragent Traction (MPa)",fontsize=20)
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                    plt.savefig(save_folder_output_png_file_path + "/tang_stress_strength_"+str(time_ind)+".png",bbox_inches = 'tight')
                #normal dir
                if plot_var_name_i == "normal_jump_rate":
                    plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all,'b-',linewidth=3.0)
                    plt.grid(False)
                    if fault_name == "mf180":
                        plt.xlim([-30000,50000])
                        plt.ylim([-0.1,0.1])
                    elif fault_name == "mf30":
                        plt.xlim([-95532.164, -30000])
                        plt.ylim([-1.0,1.0])
                    else:
                        plt.xlim([-33921.924,-160.033])
                        plt.ylim([-0.1,0.1])
                    plt.xlabel("x coord (m)",fontsize=20)
                    plt.ylabel("Normal Jump Rate (m/s)",fontsize=20)
                    if fault_name == "mf180":
                        plt.title("Main Fault Normal Jump Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Normal Jump Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Normal Jump Rate at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/normaljumprate_"+str(time_ind)+".png",bbox_inches = 'tight')
                elif plot_var_name_i == "normal_jump":
                    plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all,'b-',linewidth=3.0)
                    plt.grid(False)
                    if fault_name == "mf180":
                        plt.xlim([-30000,50000])
                        plt.ylim([-0.25,0.25])
                    elif fault_name == "mf30":
                        plt.xlim([-95532.164, -30000])
                        plt.ylim([-10.0,10.0])
                    else:
                        plt.xlim([-33921.924,-160.033])
                        plt.ylim([-0.1,0.1])
                    plt.xlabel("x coord (m)",fontsize=20)
                    plt.ylabel("Jump (m)",fontsize=20)
                    if fault_name == "mf180":
                        plt.title("Main Fault Normal Jump at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Normal Jump at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Normal Jump at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/normaljump_"+str(time_ind)+".png",bbox_inches = 'tight')
                elif plot_var_name_i == "normal_traction":
                    plt.plot(sorted_arr_czm_elem_arclength_i_all,-sorted_arr_czm_param_i_all/1e6 + dict_ini_normal_sts[fault_name],'b-',linewidth=3.0)
                    plt.grid(False)
                    if fault_name == "mf180":
                        plt.xlim([-30000,50000])
                        plt.ylim([dict_ini_normal_sts[fault_name]-500,dict_ini_normal_sts[fault_name]+500])
                    elif fault_name == "mf30":
                        plt.xlim([-95532.164, -30000])
                        plt.ylim([dict_ini_normal_sts[fault_name]-100,dict_ini_normal_sts[fault_name]+100])
                    else:
                        plt.xlim([-33921.924,-160.033])
                        plt.ylim([dict_ini_normal_sts[fault_name]-25,dict_ini_normal_sts[fault_name]+25])
                    plt.xlabel("x coord (m)",fontsize=20)
                    plt.ylabel("Normal Traction (MPa)",fontsize=20)
                    if fault_name == "mf180":
                        plt.title("Main Fault Normal Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    elif fault_name == "mf30":
                        plt.title("Fault Beyond Kink Normal Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    else:
                        plt.title("Splay Fault Normal Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.savefig(save_folder_output_png_file_path + "/normaltraction_"+str(time_ind)+".png",bbox_inches = 'tight')
                elif plot_var_name_i == "traction_y" and additional_flags == None:
                    if fault_name == "mf180":
                        plt.plot(sorted_arr_czm_elem_arclength_i_all,sorted_arr_czm_param_i_all/1e6 + dict_ini_normal_sts[fault_name] ,'b*-',linewidth=3.0)
                        plt.ylim([dict_ini_normal_sts[fault_name]-40,dict_ini_normal_sts[fault_name]+40])
                    if fault_name == "mf180":
                        plt.title("Main Fault Normal Traction at Time Step "+str(time_i)+" s",fontsize=20) #name needed to be generalized
                    plt.xlabel("Arc Length (m)",fontsize=20)
                    plt.ylabel("Normal Traction (MPa)",fontsize=20)
                    plt.xticks(fontsize=15)
                    plt.yticks(fontsize=15)
                    plt.savefig(save_folder_output_png_file_path + "/normaltraction_"+str(time_ind)+".png",bbox_inches = 'tight')
                #plt.show()
                #front position over time
                #[[[IMPORTANT: here the elements are all within the same block id, things are different if complex network fault involved]]]
                if (plot_var_name_i == "tangent_jump_rate" and additional_flags == "front_position") or (plot_var_name_i == "jump_x_rate" and additional_flags == "front_position"):
                    
                   #loop x coord of node bounds
                   #can also do y coord
                   for node_pair_index in range(len(node_pair_list_x)):
                       
                       #get current two node x coord
                       node_pair_x = node_pair_list_x[node_pair_index]
                       node_pair_y = node_pair_list_y[node_pair_index]


                       #get node_i node_j (x coord) #from node_i -> node_j
                       node_i_x = node_pair_x[0]
                       node_j_x = node_pair_x[1]
                       node_i_y = node_pair_y[0]
                       node_j_y = node_pair_y[1]

                       #arc length 
                       #only works for bf30
                       if fault_name == "bf30":
                            node_i = -1 * np.sqrt(node_i_x ** 2 + node_i_y ** 2)
                            node_j = -1 * np.sqrt(node_j_x ** 2 + node_j_y ** 2)
                       elif fault_name == "mf180":
                            node_i = node_i_x
                            node_j = node_j_x
                        
                       #find which x_i ? x_j
                       if node_i > node_j:
                           max_arclencoord = node_i
                           min_arclencoord = node_j
                           dir_flag = "decrease"
                       else:
                           max_arclencoord = node_j
                           min_arclencoord = node_i
                           dir_flag = "increase"

                       #find index of available elements within the range
                       range_index = np.where(np.logical_and(sorted_arr_czm_elem_arclength_i_all >= min_arclencoord, sorted_arr_czm_elem_arclength_i_all <= max_arclencoord))[0]

                       #find cooresponding value within the range
                       coord_range = sorted_arr_czm_elem_arclength_i_all[range_index]
                       data_range = sorted_arr_czm_param_i_all[range_index]

                       #threshold for rupture tip slip rate (m/s)
                       #data_threshold = 0.01
                       data_threshold  = 0.01

                       #find index of available elements greater than data_threshold
                       path_index = np.where(data_range > data_threshold)[0]

                       #find cooresponding x coord that is most closet to the end node_j
                       if len(path_index) != 0:
                           coord_path = coord_range[path_index]
                           front_position = closest(coord_path,node_j)
                       else:
                           front_position = node_i
                    
                       #get previous 
                       prefront = node_i
                       if node_pair_index == 0 and time_i > 0.0:
                          prefront = list_front_xcoord0[-1]
                       elif node_pair_index == 1 and time_i > 0.0:
                          prefront = list_front_xcoord1[-1]
                       
                       #if going back, make it the previous value
                       if dir_flag == "decrease" and front_position > prefront:
                            front_position = prefront
                       
                       if dir_flag == "increase" and front_position < prefront:
                            front_position = prefront
                       
                        #    #test
                        #    if time_i != 0.0:
                        #     print("coord_path:", coord_path)
                        #     print("front_position:",front_position)

                       #save into txt file with time
                       with open(save_folder_output_file_path + "/front_time"+str(node_pair_index)+".txt", "ab") as f:
                            np.savetxt(f, [time_i])
                       with open(save_folder_output_file_path + "/front_xcoord"+str(node_pair_index)+".txt", "ab") as f:
                            np.savetxt(f, [front_position])
                       
                       #save into list
                       if node_pair_index == 0:
                           list_front_time0.append(time_i)
                           list_front_xcoord0.append(front_position)
                       elif node_pair_index == 1:
                           list_front_time1.append(time_i)
                           list_front_xcoord1.append(front_position)
                # plt.close()

            # plt.savefig(save_folder_output_png_file_path + "/slip_"+str(time_ind)+".png")
            plt.close()

        if plot_var_name_i == "tangent_jump":
        
            numofrows = int(len(list_slip)/2)
            arr_slip = np.array(list_slip)
            arr_slip = arr_slip.reshape((numofrows,2))

            np.savetxt(save_folder_output_png_file_path+"/avgmaxslip"+fault_name+".txt",arr_slip,fmt='%.5f',newline=" ")

def DecodeName(exodus_file_path, decodeflag,save_folder_output_file_path):

    """
    Get Element Var Name 
    """

    #Read File
    nc = netCDF4.Dataset(exodus_file_path)

    #Get numpy.bytes name array
    if decodeflag == "name_elem_var":
        arr_name_elem_var = nc.variables['name_elem_var']
    elif decodeflag == "eb_names":
        arr_name_elem_var = nc.variables['eb_names']

    #Get number of name
    num_names = np.shape(arr_name_elem_var)[0]

    #initialize list of names
    list_name_elem_var = []

    #loop over each name
    for name_ind in range(num_names):

        name_i_decode = ''

        name_i = arr_name_elem_var[name_ind]

        # print(name_i)

        for ind in range(len(name_i)):
            
            if not np.ma.is_masked(name_i[ind]):
                
                byte_i_decode = name_i[ind].decode('UTF-8')

                name_i_decode += byte_i_decode

        list_name_elem_var.append(name_i_decode)
    
    #save
    if decodeflag == "name_elem_var":
        np.savetxt(save_folder_output_file_path + "/list_name_elem_var.txt",list_name_elem_var,fmt='%s',newline=" ")
    elif decodeflag == "eb_names":
        np.savetxt(save_folder_output_file_path + "/list_eb_names.txt",list_name_elem_var,fmt='%s',newline=" ")

def closest(lst, K):
     
    lst = lst.tolist()
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]

def CalcPlasticWork(save_folder_output_file_path):

    ##
    #Initialize dict storing elem var name / index pair
    dict_evn_index = {}

    #Read Elem Var Name
    evn = np.loadtxt(save_folder_output_file_path + "/list_name_elem_var.txt",dtype=str)

    #fill in dict_evn_index
    for evn_ind in range(len(evn)):
        evn_i = evn[evn_ind]
        # start from 1
        dict_evn_index[evn_i] = evn_ind + 1

    ##
    #Read Exodus File
    nc = netCDF4.Dataset(exodus_file_path)

    #Obtain Time Series
    timeseries = nc.variables['time_whole']

    plot_snapshot = False

    block_num = 2 #1 or 2
    B_in_data = nc.variables['vals_elem_var'+str(dict_evn_index["B_in"])+'eb'+str(block_num)][:]
    sts_change_11_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_change_11_aux"])+'eb'+str(block_num)][:]
    sts_change_12_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_change_12_aux"])+'eb'+str(block_num)][:] 
    sts_change_22_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_change_22_aux"])+'eb'+str(block_num)][:]
    sts_change_33_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_change_33_aux"])+'eb'+str(block_num)][:]
    sts_initial_11_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_initial_11_aux"])+'eb'+str(block_num)][:] 
    sts_initial_12_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_initial_12_aux"])+'eb'+str(block_num)][:]
    sts_initial_22_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_initial_22_aux"])+'eb'+str(block_num)][:]
    sts_initial_33_aux_data = nc.variables['vals_elem_var'+str(dict_evn_index["sts_initial_33_aux"])+'eb'+str(block_num)][:]

    #get t = 1.4s
    if plot_snapshot == True:
        time_slice = 14
        B_in_data = B_in_data[time_slice,:]
        sts_change_11_aux_data = sts_change_11_aux_data[time_slice,:]
        sts_change_12_aux_data = sts_change_12_aux_data[time_slice,:]
        sts_change_22_aux_data = sts_change_22_aux_data[time_slice,:]
        sts_change_33_aux_data = sts_change_33_aux_data[time_slice,:]
        sts_initial_11_aux_data = sts_initial_11_aux_data[time_slice,:]
        sts_initial_12_aux_data = sts_initial_12_aux_data[time_slice,:]
        sts_initial_22_aux_data = sts_initial_22_aux_data[time_slice,:]
        sts_initial_33_aux_data = sts_initial_33_aux_data[time_slice,:]

    #compute total stress components
    sts_total_11_aux = sts_change_11_aux_data + sts_initial_11_aux_data
    sts_total_12_aux = sts_change_12_aux_data + sts_initial_12_aux_data
    sts_total_22_aux = sts_change_22_aux_data + sts_initial_22_aux_data
    sts_total_33_aux = sts_change_33_aux_data + sts_initial_33_aux_data

    #compute trace of total stress
    sts_total_trace = sts_total_11_aux + sts_total_22_aux + sts_total_33_aux

    #compute deviatroic stress components
    sts_d_11 = sts_total_11_aux - 1.0/3.0 * sts_total_trace
    sts_d_12 = sts_total_12_aux
    sts_d_22 = sts_total_22_aux - 1.0/3.0 * sts_total_trace
    sts_d_33 = sts_total_33_aux - 1.0/3.0 * sts_total_trace

    #compute plastic strain rate
    Cg = 1e-10
    m1 = 10
    eps_p_rate_11 = Cg * B_in_data ** m1 * sts_d_11
    eps_p_rate_12 = Cg * B_in_data ** m1 * sts_d_12
    eps_p_rate_22 = Cg * B_in_data ** m1 * sts_d_22
    eps_p_rate_33 = Cg * B_in_data ** m1 * sts_d_33

    print(np.max(eps_p_rate_11))

    #compute plastic work rate
    plastic_work_rate = eps_p_rate_11 * sts_total_11_aux + 2 * eps_p_rate_12 * sts_total_12_aux + eps_p_rate_22 * sts_total_22_aux + eps_p_rate_33 * sts_total_33_aux

    print(np.shape(plastic_work_rate))

    if plot_snapshot == False:
        plastic_work_rate_t_list = np.sum(plastic_work_rate, axis=1)
        plastic_work_rate_append_list = []
        for i in range(len(plastic_work_rate_t_list)):
            plastic_work_rate_append_list.append(plastic_work_rate_t_list[i] + sum(plastic_work_rate_t_list[:i]))
        for j in range(len(plastic_work_rate_append_list)):
            plastic_work_rate_append_list[j] = plastic_work_rate_append_list[j] * 0.01 * 2 * 540

    #save res
    np.savetxt(save_folder_output_file_path+"/plastic_work/timeseries.txt",timeseries,fmt='%.5f',newline=" ")
    np.savetxt(save_folder_output_file_path+"/plastic_work/plastic_work.txt",plastic_work_rate_append_list,fmt='%.5f',newline=" ")

    plt.figure()
    plt.plot(plastic_work_rate_append_list)
    plt.show()

    exit()

if __name__ == "__main__":

    #file path (read)
    mf180_ptrs_data_file_path = "./ptrsdata/mf180_ptrs_data.txt"

    exodus_file_path = "/Users/andyz/projects/farms/examples/cdbm_planarfault_long3/Cd1e7local/restart_planarfault_main_out.e"

    #file path (save)
    mf180_save_folder_output_file_path = "./outputs"

    mf180_save_folder_output_png_file_path = "./outputs"

    ##
    dict_save_folder_output_file_path = {"mf180" : mf180_save_folder_output_file_path}
    
    dict_save_folder_output_png_file_path = {"mf180" : mf180_save_folder_output_png_file_path}

    ##
    dict_ptrs_data_file_path = { "mf180" : mf180_ptrs_data_file_path}

    ##
    #line plot 2D
    # plot_var_name = ['tangent_jump_rate', 'tangent_traction','tangent_jump']
    # plot_var_name = ['tangent_traction']
    # plot_var_name = ["normal_jump_rate", "normal_jump", "normal_traction"]
    # plot_var_name = ["tangent_jump_rate"]
    # plot_var_name = ["jump_x"]
    # plot_var_name = ["traction_x"]
    # plot_var_name = ["traction_y"]
    plot_var_name = ["tangent_jump_rate", "tangent_jump"]
    additional_files = {"mf180" : [(0,0),(10000,0),(0,0),(-10000,0)]}
    
    # additional_flags = "front_position"
    # additional_flags = "stress_and_strength"
    additional_flags = None

    nc = netCDF4.Dataset(exodus_file_path)
    # print(nc.variables['name_elem_var'])
    # print(np.shape(nc.variables['time_whole']))
    # print(len(nc.variables['eb_status']))
    # print(nc.variables['eb_names'])

    ##
    list_fault_name = ["mf180"]

    #retrieve element/node info for CZM
    #prior knowledge of master locs
    master_locs_mf180 = [0,1,0,0]              # "0" * x + "1" * y + "0" ">" 0

    dict_master_locs = {"mf180" : master_locs_mf180}

    run_decode_retrieve_flag = False

    if run_decode_retrieve_flag == True:
        for fault_name_index in range(len(list_fault_name)):

            #current fault name
            fault_name_i = list_fault_name[fault_name_index]

            print("Current Fault: ",fault_name_i)

            #get element names
            DecodeName(exodus_file_path,"name_elem_var",dict_save_folder_output_file_path[fault_name_i])

            #get element block names
            DecodeName(exodus_file_path,"eb_names",dict_save_folder_output_file_path[fault_name_i])

            #retrieve element/node info for CZM
            RetrieveCZMElementIDs(exodus_file_path = exodus_file_path,
                                save_list_czm_elem_file_path = dict_save_folder_output_file_path[fault_name_i],
                                dict_ptrs_data_file_path = dict_ptrs_data_file_path[fault_name_i],
                                master_locs = dict_master_locs[fault_name_i],
                                fault_name=fault_name_i)

    #line plot 2D
    fault_name = "mf180"

    #initial shear stress
    dict_ini_shear_sts = {"mf180" : 70}
    
    #initial normal stress
    dict_ini_normal_sts = {"mf180" : 120}

    #plastic work
    CalcPlasticWork(save_folder_output_file_path=dict_save_folder_output_file_path[fault_name])

    exit()
    LinePlotParamVal2D(exodus_file_path=exodus_file_path,
                        plot_var_name=plot_var_name,
                        save_folder_output_file_path=dict_save_folder_output_file_path[fault_name],
                        save_folder_output_png_file_path=dict_save_folder_output_png_file_path[fault_name],
                        fault_name=fault_name,
                        additional_files=additional_files,
                        additional_flags=None,
                        dict_ini_shear_sts=dict_ini_shear_sts,
                        dict_ini_normal_sts=dict_ini_normal_sts)


    ##
    #calculate front position vs time
    additional_flags = "front_position"
    plotfront = True
    calcfront = True

    if calcfront == True:
        LinePlotParamVal2D(exodus_file_path=exodus_file_path,
                        plot_var_name=plot_var_name,
                        save_folder_output_file_path=dict_save_folder_output_file_path[fault_name],
                        save_folder_output_png_file_path=dict_save_folder_output_png_file_path[fault_name],
                        fault_name=fault_name,
                        additional_files=additional_files,
                        additional_flags=additional_flags)
    if plotfront == True:
        # print front position vs time ##
        # 1
        datatime1s05 = np.loadtxt(dict_save_folder_output_file_path[fault_name] + "/front_time0.txt")
        datapos1s05 = np.loadtxt(dict_save_folder_output_file_path[fault_name] + "/front_xcoord0.txt")
        # 2
        datatime1s052 = np.loadtxt(dict_save_folder_output_file_path[fault_name] + "/front_time1.txt")
        datapos1s052 = np.loadtxt(dict_save_folder_output_file_path[fault_name] + "/front_xcoord1.txt")
        plt.figure()
        plt.plot(datatime1s05,datapos1s05,'r-')
        plt.plot(datatime1s052,datapos1s052,'r-')
        if fault_name == "bf30":
            plt.title("splay fault front position vs time (junction point dir)")
            plt.xlim([0,20])
        elif fault_name == "mf180":
            plt.title("main fault front position vs time (right rupture dir)")
        plt.xlabel("time")
        plt.ylabel("front position")
        plt.legend()
        plt.show()