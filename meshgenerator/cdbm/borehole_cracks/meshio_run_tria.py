#For TRIA3
import matplotlib.pyplot as plt
import numpy as np
import math
from labellines import labelLines
import gmsh
import pandas as pd
import meshio
import csv
import os
from tqdm import tqdm

m = meshio.read('./borehole_wofaults.msh')

def find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y):

    #fault line slope 
    slope = (end_ptr_y - start_ptr_y) / (end_ptr_x - start_ptr_x)

    #prep line slope
    if slope != 0:
        slope_prep = - 1 / slope
    else:
        slope_prep = 0

    #segement 
    b = start_ptr_y - slope * start_ptr_x

    #segement (prep + start_ptr)
    b_start_prep = start_ptr_y - slope_prep * start_ptr_x

    #segement (prep + end_ptr)
    b_end_prep = end_ptr_y - slope_prep * end_ptr_x

    return slope,slope_prep,b,b_start_prep, b_end_prep

#---------------------------#
# Physical Lines            #
#---------------------------#
#number of lines
numoflines = 12

phylines = np.zeros((numoflines,4))

phylines[0,:] = [200, 0, 800, 0]; 
phylines[1,:] = [173.20508, 100, 692.82032, 400]; 
phylines[2,:] = [100, 173.20508, 400, 692.82032];
phylines[3,:] = [0, 200, 0, 800]; 
phylines[4,:] = [-100, 173.20508, -400, 692.82032]; 
phylines[5,:] = [-173.20508, 100, -692.82032, 400];
phylines[6,:] = [-200, 0, -800, 0];
phylines[7,:] = [-173.20508, -100, -692.82032, -400];
phylines[8,:] = [-100, -173.20508, -400, -692.82032];
phylines[9,:] = [0, -200, 0, -800];
phylines[10,:] = [100, -173.20508, 400, -692.82032];
phylines[11,:] = [173.20508, -100, 692.82032, -400];

#generate lists
try:
    os.remove('./outputs/list_elemid_embeded.txt')
    os.remove('./outputs/list_subdomainid_embeded.txt')
except OSError:
    pass
with open('./outputs/list_elemid_embeded.txt', 'w'):
    pass
with open('./outputs/list_subdomainid_embeded.txt', 'w'):
    pass

for i in tqdm(range(1,np.shape(phylines)[0]+1)): #assume embedded line name from 1

    #get coordinates
    ptr_i_x = phylines[i-1,0]
    ptr_i_y = phylines[i-1,1]
    ptr_j_x = phylines[i-1,2]
    ptr_j_y = phylines[i-1,3]

    #not align with x or y coordinate
    if (ptr_j_y - ptr_i_y != 0.0 and ptr_j_x - ptr_i_x != 0.0 ):
        
        #fault line slope 
        slope = (ptr_j_y - ptr_i_y) / (ptr_j_x - ptr_i_x)
        
        # determine start/end ptr based on slope
        if slope < 0:
            if ( ptr_i_x < ptr_j_x ): #min x ptr is start_ptr
                start_ptr_x = ptr_i_x
                start_ptr_y = ptr_i_y
                end_ptr_x = ptr_j_x
                end_ptr_y = ptr_j_y
            else:
                start_ptr_x = ptr_j_x
                start_ptr_y = ptr_j_y
                end_ptr_x = ptr_i_x
                end_ptr_y = ptr_i_y
        elif slope > 0:
            if ( ptr_i_x < ptr_j_x ): #max x ptr is start_ptr
                start_ptr_x = ptr_j_x
                start_ptr_y = ptr_j_y
                end_ptr_x = ptr_i_x
                end_ptr_y = ptr_i_y
            else:
                start_ptr_x = ptr_i_x
                start_ptr_y = ptr_i_y
                end_ptr_x = ptr_j_x
                end_ptr_y = ptr_j_y
        
        #find bounds on lines
        slope,slope_prep,b,b_start,b_end = find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y)

        # Get Physical Line Element Indices
        phy_nodal_ind = m.cell_sets_dict['embeded'+str(i)]['line']
        #print(phy_nodal_ind)

        # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
        phy_nodal_arr = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind,:]))
        #print(phy_nodal_arr)

        # Get QUAD4 Element Connectivity
        tria_elem_connect = m.cells_dict['triangle']

        # Find Elements Along the Crack Line
        ind_elem_crack_upper = []
        ind_elem_crack_lower = []
        num_elem = np.shape(tria_elem_connect)[0] # elem num
        for elem_ind in range(num_elem):   
            elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
            num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr,elem_connect_i)) # count
            if (num_ptr_collapse >= 1): # this is a crack aligned element
                coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
                coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
                if ( len(np.argwhere(slope_prep * coord_data_x + b_start - coord_data_y >= 0 )) > 1 and len(np.argwhere(slope_prep * coord_data_x + b_end - coord_data_y <= 0)) > 1 ): #slope_prep > 0
                    if ( slope * np.sum(coord_data_x)/3 + b - np.sum(coord_data_y)/3 < 0 ): #hardcode #above the fault #center ptr
                        ind_elem_crack_upper.append(elem_ind)
                    else:
                        ind_elem_crack_lower.append(elem_ind)
                else:
                    continue

        #only store branch
        ind_elem_all_embeded = ind_elem_crack_upper + ind_elem_crack_lower

        ind_subdomainid_all_embeded = []
        for k in range(len(ind_elem_crack_upper)):
            ind_subdomainid_all_embeded.append(100+i)
        for l in range(len(ind_elem_crack_lower)):
            ind_subdomainid_all_embeded.append(200+i)

        #save elem ids
        with open("./outputs/list_elemid_embeded.txt", "ab") as f:
            np.savetxt(f, ind_elem_all_embeded, fmt='%i', newline=" ")

        with open("./outputs/list_subdomainid_embeded.txt", "ab") as f:
            np.savetxt(f, ind_subdomainid_all_embeded, fmt='%i', newline=" ")

    else:

        #determine the case
        if ptr_j_x - ptr_i_x == 0.0: #along ycoord
            
            if ptr_i_y < ptr_j_y: 
                start_ptr = ptr_i_y
                end_ptr = ptr_j_y
            else:
                start_ptr = ptr_j_y
                end_ptr = ptr_i_y               
        
        else: #along xcoord
            
            if ptr_i_x < ptr_j_x:
                start_ptr = ptr_i_x
                end_ptr = ptr_j_x  
            else:
                start_ptr = ptr_j_x
                end_ptr = ptr_i_x                        

        # Get Physical Line Element Indices
        phy_nodal_ind = m.cell_sets_dict['embeded'+str(i)]['line']
        #print(phy_nodal_ind)

        # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
        phy_nodal_arr = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind,:]))
        #print(phy_nodal_arr)

        # Get QUAD4 Element Connectivity #Need to define Physical Surface
        tria_elem_connect = m.cells_dict['triangle']

        # Find Elements Along the Crack Line
        ind_elem_crack_upper = []
        ind_elem_crack_lower = []
        num_elem = np.shape(tria_elem_connect)[0] # elem num
        for elem_ind in range(num_elem):   
            elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
            #print(elem_connect_i)
            num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr,elem_connect_i)) # count
            if (num_ptr_collapse >= 1): # this is a crack aligned element
                coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
                coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
                
                if ptr_j_x - ptr_i_x == 0.0: #along y coordinate
                    coord_data_para = coord_data_y
                    coord_data_prep = coord_data_x 
                else: #along x coordinate
                    coord_data_para = coord_data_x
                    coord_data_prep = coord_data_y    

                if ( len(np.argwhere(coord_data_para >= start_ptr)) > 1 and len(np.argwhere(coord_data_para < end_ptr)) > 1 ): #hardcode #within fault region
                    if ( np.sum(coord_data_prep) > 0 ): #hardcode #above the fault
                        ind_elem_crack_upper.append(elem_ind)
                    else:
                        ind_elem_crack_lower.append(elem_ind)
                else:
                    continue
        
        #only store branch
        ind_elem_all_embeded = ind_elem_crack_upper + ind_elem_crack_lower

        ind_subdomainid_all_embeded = []
        for k in range(len(ind_elem_crack_upper)):
            ind_subdomainid_all_embeded.append(100+i)
        for l in range(len(ind_elem_crack_lower)):
            ind_subdomainid_all_embeded.append(200+i)
        
        #save elem ids
        with open("./outputs/list_elemid_embeded.txt", "ab") as f:
            np.savetxt(f, ind_elem_all_embeded, fmt='%i', newline=" ")

        with open("./outputs/list_subdomainid_embeded.txt", "ab") as f:
            np.savetxt(f, ind_subdomainid_all_embeded, fmt='%i', newline=" ")
