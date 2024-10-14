import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import os
import glob
from tqdm import tqdm
from scipy.spatial import distance
import shutil

class postprocessclass:

    def __init__(self, 
                 file_path = None, 
                 decodeflags = None, 
                 plotvar = None, 
                 plotvar_nodal = None,
                 save_file_path = None,
                 dim = None,
                 node_per_elem = None):

        """
        class initialization
        output: self.file_path 
                self.decodeflags 
                self.plotvar 
                self.save_file_path
                self.dim
                self.node_per_elem

        other quantities: self.
        """

        #read exodus file
        self.nc = netCDF4.Dataset(file_path)
        self.decodeflags = decodeflags
        self.plotvar = plotvar
        self.plotvar_nodal = plotvar_nodal
        self.save_file_path = save_file_path
        self.dim = dim
        self.node_per_elem = node_per_elem

    def decode_name(self):

        """
        Get Element/Node Var Name 
        output: self.evn
                self.dict_evn_index (elem)
                self.dict_evn_index_nodal (node)
        """
        
        #Initialize dict storing elem var name / index pair
        self.dict_evn_index = {}

        #Initialize dict storing elem var name / index pair
        self.dict_evn_index_nodal = {}

        #loop over all decodeflags
        for i in range(len(self.decodeflags)):

            decodeflag_i = self.decodeflags[i]

            print("Decode "+decodeflag_i+" ...")

            #Get numpy.bytes name array
            if decodeflag_i == "name_elem_var":
                arr_name_elem_var = self.nc.variables['name_elem_var']
            elif decodeflag_i == "name_nod_var":
                arr_name_elem_var = self.nc.variables['name_nod_var']
            elif decodeflag_i == "eb_names":
                arr_name_elem_var = self.nc.variables['eb_names']

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
            if decodeflag_i == "name_elem_var":
                np.savetxt(self.save_file_path + "/list_name_elem_var.txt",list_name_elem_var,fmt='%s',newline=" ")
                #load the names and return 
                self.evn = np.loadtxt(self.save_file_path + "/list_name_elem_var.txt",dtype=str)
                
                #fill in dict_evn_index
                for evn_ind in range(len(self.evn)):
                    evn_i = self.evn[evn_ind]
                    # start from 1
                    self.dict_evn_index[evn_i] = evn_ind + 1

            elif decodeflag_i == "name_nod_var":
                np.savetxt(self.save_file_path + "/list_name_nod_var.txt",list_name_elem_var,fmt='%s',newline=" ")
                #load the names and return 
                self.evn_nodal = np.loadtxt(self.save_file_path + "/list_name_nod_var.txt",dtype=str)
                
                if np.shape(self.evn_nodal)!=(): #multiple variables

                    #fill in dict_evn_index
                    for evn_ind in range(len(self.evn_nodal)):
                        evn_i = self.evn_nodal[evn_ind]
                        # start from 1
                        self.dict_evn_index_nodal[evn_i] = evn_ind + 1

                else: #there is only one variable
        
                    evn_i = self.evn_nodal
                    self.dict_evn_index_nodal[str(evn_i)] = 1
            

            elif decodeflag_i == "eb_names":
                np.savetxt(self.save_file_path + "/list_eb_names.txt",list_name_elem_var,fmt='%s',newline=" ")
                #load the names and return 
                self.evn = np.loadtxt(self.save_file_path + "/list_name_elem_var.txt",dtype=str)

    def pre_post(self,
                 elemblock_ind = None):

        """
        Get element info
        output: self.var_ind
                self.num_elem_block
                self.list_elem_dists
                self.num_elem
        """

        #get total number of block inds
        self.num_elem_block = len(self.nc.variables['eb_status'])

        #if specific block is given, then no need to loop over all blocks
        if elemblock_ind is not None:

            #get connectivity for the current block
            #shrift by 1
            elem_connect = self.nc.variables['connect' + str(elemblock_ind)][:]-1

            #TET4 Elem Num
            self.num_elem = np.shape(elem_connect)[0]

            #Initialize list for saving centroid point of Tet4 element
            self.list_elem_centroids = np.zeros((self.num_elem,self.dim))

            #Loop over elements
            print("Get Element Info ...")
            for elem_ind in tqdm(range(self.num_elem)):

                #get connectivity of current element
                elem_connect_i = elem_connect[elem_ind,:]  

                #get coordinate for current element
                coord_data_x = self.nc.variables['coordx'][elem_connect_i]
                coord_data_y = self.nc.variables['coordy'][elem_connect_i]  
                if self.dim == 3:
                    coord_data_z = self.nc.variables['coordz'][elem_connect_i] 

                #get centroid by taking average of each dimension
                self.list_elem_centroids[elem_ind,0] = np.sum(coord_data_x)/self.node_per_elem
                self.list_elem_centroids[elem_ind,1] = np.sum(coord_data_y)/self.node_per_elem
                if self.dim == 3:
                    self.list_elem_centroids[elem_ind,2] = np.sum(coord_data_z)/self.node_per_elem

                

        else:

            #loop over block ind
            for elemblock_ind in range(1,self.num_elem_block+1):

                #get connectivity for the current block
                #shrift by 1
                elem_connect = self.nc.variables['connect' + str(elemblock_ind)][:]-1

                #TET4 Elem Num
                self.num_elem = np.shape(elem_connect)[0]

                #Initialize list for saving centroid point of Tet4 element
                self.list_elem_centroids = np.zeros((self.num_elem,self.dim))

                #Loop over elements
                print("Get Element Info ...")
                for elem_ind in tqdm(range(self.num_elem)):

                    #get connectivity of current element
                    elem_connect_i = elem_connect[elem_ind,:]  

                    #get coordinate for current element
                    coord_data_x = self.nc.variables['coordx'][elem_connect_i]
                    coord_data_y = self.nc.variables['coordy'][elem_connect_i]  
                    if self.dim == 3:
                        coord_data_z = self.nc.variables['coordz'][elem_connect_i] 

                    #get centroid by taking average of each dimension
                    self.list_elem_centroids[elem_ind,0] = np.sum(coord_data_x)/self.node_per_elem
                    self.list_elem_centroids[elem_ind,1] = np.sum(coord_data_y)/self.node_per_elem
                    if self.dim == 3:
                        self.list_elem_centroids[elem_ind,2] = np.sum(coord_data_z)/self.node_per_elem

    def post_nodal(self, 
                   ptr_coord,
                   custom_name,
                   ptr_index):

        """
        This function is used to postprocessing nodal quantities
        """

        print(self.dict_evn_index_nodal)

        #get required param to be plotted
        var_ind = self.dict_evn_index_nodal[self.plotvar_nodal]
        
        print("Loading val array ...")

        #load param array
        arr_param = self.nc.variables['vals_nod_var'+str(var_ind)][:,:]

        #current point x y z coordinate
        ptr_x = ptr_coord[0]
        ptr_y = ptr_coord[1]
        ptr_z = ptr_coord[2]

        print("Loading coordinate array ...")

        #load coordinate
        x_coord = self.nc.variables['coordx']
        y_coord = self.nc.variables['coordy']
        z_coord = self.nc.variables['coordz']

        #find index of that point
        idc = self.helper_find_nearest(x_coord[:], ptr_x, y_coord[:], ptr_y, z_coord[:], ptr_z)

        print("Loading time array ...")

        #Obtain Time Series
        timeseries = self.nc.variables['time_whole']

        print("Slicing the x,y,z values ...")

        #get time history of current ptr
        ptr_valhist = arr_param[:,idc]

        #save data
        np.savetxt(self.save_file_path+"/list_ptr_"+str(ptr_index)+"_"+self.plotvar_nodal+"_of_"+custom_name+".txt",ptr_valhist)
        np.savetxt(self.save_file_path+"/list_timeseries.txt",timeseries)

    def post_elemental(self,
                       ptr_coord,
                       elemblock_ind):
        
        """
        This function is used to postprocessing nodal quantities
        """

        for plotvar_i in range(len(self.plotvar)):

            #get required param to be plotted
            var_ind = self.dict_evn_index[self.plotvar[plotvar_i]]
            
            print("Loading val array ...")

            #load param array
            arr_param = self.nc.variables['vals_elem_var'+str(var_ind)+'eb'+str(elemblock_ind)][:,:]

            #current point x y z coordinate
            ptr_x = ptr_coord[0]
            ptr_y = ptr_coord[1]
            ptr_z = ptr_coord[2]

            print("Loading coordinate array ...")

            #load coordinate
            x_coord = self.list_elem_centroids[:,0]
            y_coord = self.list_elem_centroids[:,1]
            if self.dim == 3:
                z_coord = self.list_elem_centroids[:,2]

            #find index of that point
            if self.dim == 3:
                idc = self.helper_find_nearest(x_coord[:], ptr_x, y_coord[:], ptr_y, z_coord[:], ptr_z)
            else:
                idc = self.helper_find_nearest(x_coord[:], ptr_x, y_coord[:], ptr_y)

            print("Loading time array ...")

            #Obtain Time Series
            timeseries = self.nc.variables['time_whole']

            print("Slicing the x,y,z values ...")

            #get time history of current ptr
            ptr_valhist = arr_param[:,idc]

            #save data
            if self.dim == 3:
                np.savetxt(self.save_file_path+'/'+self.plotvar[plotvar_i]+'_strike'+str(ptr_x/1000)+'_dip'+str(ptr_y/1000)+'.txt',
                    ptr_valhist,
                    fmt='%.7f',
                    newline=" ")
            else:
                np.savetxt(self.save_file_path+'/'+self.plotvar[plotvar_i]+'_strike'+str(ptr_x/1000)+'_dip'+str(ptr_z/1000)+'.txt',
                    ptr_valhist,
                    fmt='%.7f',
                    newline=" ")               
            np.savetxt(self.save_file_path+"/list_timeseries.txt",timeseries)
            np.savetxt(self.save_file_path+"/list_elem_centroids.txt",self.list_elem_centroids)

    def helper_find_nearest(self,array_x, value_x, array_y, value_y, array_z = None, value_z = None):
    
        """
        This function is used to find the closet node to the given coordinates
        """

        print("Start finding nearest point ...")

        if self.dim == 3:
            
            #
            array_x = array_x[:].reshape((len(array_x),1))
            array_y = array_y[:].reshape((len(array_y),1))
            array_z = array_z[:].reshape((len(array_y),1))

            #
            nodes = np.hstack((array_x,array_y,array_z))
            node = np.hstack((value_x,value_y,value_z))

            closest_index = distance.cdist([node], nodes).argmin()

        else:

            #
            array_x = array_x[:].reshape((len(array_x),1))
            array_y = array_y[:].reshape((len(array_y),1))

            #
            nodes = np.hstack((array_x,array_y))
            node = np.hstack((value_x,value_y))

            closest_index = distance.cdist([node], nodes).argmin()


            print("Find point! x: "+str(nodes[closest_index][0])+" y: "+str(nodes[closest_index][1]))

        return closest_index

    ###The functions below are for plotting, should not general and are subject to change###

    def plot_damage_radially_distribution(self,
                                          custom_names,
                                          custom_styles):

        """
        Plot damage figures
        """

        #initialize figure
        plt.figure()

        #loop over cases
        for i in range(len(custom_names)):

            list_dists_interval_i = np.loadtxt(self.save_file_path+"/list_dists_interval_"+custom_names[i]+".txt")
            lists_damage_interval_percent_i = np.loadtxt(self.save_file_path+"/list_damage_interval_percent_"+custom_names[i]+".txt")

            plt.plot(list_dists_interval_i,
                     lists_damage_interval_percent_i,
                     custom_styles[i],label=custom_names[i])
        
        plt.legend()
        plt.xlabel("Relative distance from the center r (m)")
        plt.ylabel("Percent crack damage d (%)")
        plt.title("Crack Damage Radially Distribution")
        plt.savefig(self.save_file_path+"/damage_plots.png")

    def plot_pressure_timehist(self,
                                custom_names,
                                custom_styles,
                                plotvars_nodal,
                                ptr_num):
        
        """
        This function plot point-wise comparison
        """
        
        #loop over points
        for i in range(ptr_num):

            ptr_index = i

            #initialize figure
            plt.figure()

            #load time
            timeseries = np.loadtxt(self.save_file_path+"/list_timeseries.txt")

            #loop over cases
            for j in range(len(custom_names)):

                #get custom_name
                custom_name = custom_names[j]

                #get plotvar_nodal
                plotvar_nodal = plotvars_nodal[j]

                #get custom_style
                custom_style = custom_styles[j]

                #load file    
                ptr_var = np.loadtxt(self.save_file_path+"/list_ptr_"+str(ptr_index)+"_"+plotvar_nodal+"_of_"+custom_name+".txt")

                #plot
                plt.plot(timeseries*1e6,ptr_var/1e6,custom_style,label=custom_name)

            #add labels
            plt.legend()
            plt.xlabel("Time (ms)")
            plt.ylabel("Pressure (MPa)")
            plt.title("Pressure Time History")
            plt.savefig(self.save_file_path+"/pressure_plots_at_ptr"+str(i)+".png") 

    def plot_var_distribution(self,
                              plot_var_distribution,
                              scale):

        """
        This function plots the distribution of a elemental variable
        """

        #get the variable plot index
        self.var_ind = self.dict_evn_index[plot_var_distribution]

        #initialize the plot
        plt.figure()

        #loop over block ind
        for elemblock_ind in range(1,self.num_elem_block+1):

            #obtain var data array (last time step)
            arr_var = self.nc.variables['vals_elem_var'+str(self.var_ind)+'eb'+str(elemblock_ind)][-1,:] 
            np.savetxt(self.save_file_path+"/arr_var.txt",arr_var)       

            #output maximum and minimum value of the array
            print("Maximum value of variable "+plot_var_distribution+" :", np.max(arr_var)/scale)
            print("Minimum value of variable "+plot_var_distribution+" :", np.min(arr_var)/scale)

            #plot the distribution   
            plt.hist(arr_var/scale, bins=30, edgecolor='black')
            plt.title("Histogram for"+plot_var_distribution)
            plt.show()

class systemops:

    """
    This class is for some quick-os opeartions
    """

    def ensure_folder_exists(folder_path):
        # Check if the folder exists
        if not os.path.exists(folder_path):
            # If the folder does not exist, create it
            os.makedirs(folder_path)
            print(f"Created folder: {folder_path}")
        else:
            print(f"Folder already exists: {folder_path}")

    def remove_all_files_in_folder(folder_path):

        # Use glob to get all file paths in the folder
        files = glob.glob(os.path.join(folder_path, '*'))
        
        # Iterate over the file paths and remove each file
        for f in files:
            try:
                os.remove(f)
                print(f"Removed: {f}")
            except OSError as e:
                print(f"Error: {e.strerror} - {f}")

    def remove_specific_folders(folder_path, folders_to_remove):
        # Iterate over the list of folders you want to remove
        for folder in folders_to_remove:
            full_path = os.path.join(folder_path, folder)
            if os.path.isdir(full_path):  # Check if it's a directory
                try:
                    shutil.rmtree(full_path)  # Recursively delete the directory
                    print(f"Removed folder: {full_path}")
                except OSError as e:
                    print(f"Error: {e.strerror} - {full_path}")
            else:
                print(f"Folder not found: {full_path}")
    
    def rename_folder(old_name, new_name):
        try:
            os.rename(old_name, new_name)
            print(f"Folder '{old_name}' renamed to '{new_name}' successfully.")
        except OSError as e:
            print(f"Error: {e}")






            

