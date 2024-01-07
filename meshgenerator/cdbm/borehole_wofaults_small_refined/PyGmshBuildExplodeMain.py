"""
Multifault Generation

Main Executable

Created: Chunhui Zhao, Mar 29, 2023
"""

# import module 
import numpy as np
import gmsh
import math
import matplotlib.pyplot as plt
from scipy.stats import qmc
import meshio
from scipy.spatial import distance

#User-input data
#Coordinates of Corner Points
"""
case_flag: 2D-Cluster
4 - - - - - 3
|           |
|           |
|           |
1 - - - - - 2

case_flag: 2D-Cluster-Well - Hole presented in the domain as inner boundary
4 - - - - - 3
|           |
|     o     |
|           |
1 - - - - - 2
"""
#case identifier
case_flag = "2D-Cluster-Well"

# global path
global_path = "."

##file path
#file saving coarse mesh
if case_flag == "2D-Cluster":
    file_path1 = global_path + "/mshfiles/multifaults.msh"
    file_path2 = global_path + "/mshfiles/network.msh"
    txt_elems_path = global_path + "/elemdata/elem_ind.txt"
    txt_ids_path = global_path + "/elemdata/elem_ids.txt"
    txt_surround_path = global_path + "/elemdata/elem_surroundingblocks.txt"
elif case_flag == "2D-Cluster-Well":
    file_path2 = global_path + "/borehole_wofaults.msh"
    txt_elems_path = global_path + "/elemdata/elem_ind_well.txt"
    txt_ids_path = global_path + "/elemdata/elem_ids_well.txt"
    txt_surround_path = global_path + "/elemdata/elem_surroundingblocks_well.txt"
    txt_global_path = global_path + "/elemdata/"

csv_file_path = global_path + "/elemdata/area.csv"

png1_file_path = global_path + "/elemdata/1stptr.png"
png2_file_path = global_path + "/elemdata/lineplot.png"

radius_2nd = 200
lc2 = 10
    
#read file
m = meshio.read(file_path2)

#Get TRIA3 Connectivity
tria_elem_connect = m.cells_dict['triangle']

#Get coords
bound_min_x = -1000 #x1
bound_max_x = 1000 #x2
bound_min_y = -1000 #y1
bound_max_y = 1000 #y3

borehole_center_x = 0
borehole_center_y = 0

# Find Elements Along Physical Groups
## initialize containers
ind_elem_left_boundary = []
ind_elem_right_boundary = []
ind_elem_top_boundary = []
ind_elem_bottom_boundary = []
ind_elem_borehole_boundary = []
##
count = 0
num_elem = np.shape(tria_elem_connect)[0] # elem num
for elem_ind in range(num_elem):
    
    #
    count += 1
    print("find physical groups process: ", elem_ind/num_elem * 100, "%")
    
    ##   
    elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
    coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
    coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray

    #get centroid point coordinate for current element
    coord_centriod_x = np.sum(coord_data_x) / 3
    coord_centriod_y = np.sum(coord_data_y) / 3

    #add to the group
    num_ptr_collapse_leftbc   = np.count_nonzero(np.isin(coord_data_x,bound_min_x))
    num_ptr_collapse_rightbc  = np.count_nonzero(np.isin(coord_data_x,bound_max_x)) 
    num_ptr_collapse_bottombc = np.count_nonzero(np.isin(coord_data_y,bound_min_y)) 
    num_ptr_collapse_topbc    = np.count_nonzero(np.isin(coord_data_y,bound_max_y)) 
    if  ( num_ptr_collapse_leftbc   == 2 ): ind_elem_left_boundary.append(elem_ind)
    elif( num_ptr_collapse_rightbc  == 2 ): ind_elem_right_boundary.append(elem_ind)
    elif( num_ptr_collapse_bottombc == 2 ): ind_elem_bottom_boundary.append(elem_ind)
    elif( num_ptr_collapse_topbc    == 2 ): ind_elem_top_boundary.append(elem_ind)

    #compute dist and add it to group
    dist = np.sqrt( (coord_centriod_x - borehole_center_x) ** 2 + (coord_centriod_y - borehole_center_y) ** 2 )
    if ( dist < radius_2nd + lc2 and dist >= radius_2nd ): ind_elem_borehole_boundary.append(elem_ind)


# Generate subdomain IDs
ind_subdomainid_left_boundary     = [1000] * len(ind_elem_left_boundary)
ind_subdomainid_right_boundary    = [1001] * len(ind_elem_right_boundary)
ind_subdomainid_top_boundary      = [1002] * len(ind_elem_top_boundary)
ind_subdomainid_bottom_boundary   = [1003] * len(ind_elem_bottom_boundary)
ind_subdomainid_borehole_boundary = [1004] * len(ind_elem_borehole_boundary)

# Gather ids
ind_elem_all_boundary = ind_elem_left_boundary + ind_elem_right_boundary + ind_elem_top_boundary + ind_elem_bottom_boundary + ind_elem_borehole_boundary
ind_subdomainid_all_boundary = ind_subdomainid_left_boundary + ind_subdomainid_right_boundary + ind_subdomainid_top_boundary + ind_subdomainid_bottom_boundary + ind_subdomainid_borehole_boundary

#save
#save txt file
np.savetxt(txt_global_path + "elem_ind_boundaries.txt",ind_elem_all_boundary,fmt='%i',newline=" ")

np.savetxt(txt_global_path + "elem_ids_boundaries.txt",ind_subdomainid_all_boundary,fmt='%i',newline=" ")

