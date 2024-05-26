import numpy as np
import meshio

# Load Files
m = meshio.read('../TPV24.msh')

# Get TETRA Element Nodal Connectivity
tetra_elem_connect = m.cells_dict['tetra']

print(m.cell_sets)
exit()

#
def find_prep_line_eq(start_ptr_x_surf2,start_ptr_y_surf2,end_ptr_x_surf2,end_ptr_y_surf2):

    #fault line slope_surf2 
    slope_surf2 = (end_ptr_y_surf2 - start_ptr_y_surf2) / (end_ptr_x_surf2 - start_ptr_x_surf2)

    #prep line slope_surf2
    slope_surf2_prep = - 1 / slope_surf2

    #segement 
    b = start_ptr_y_surf2 - slope_surf2 * start_ptr_x_surf2

    #segement (prep + start_ptr)
    b_start_prep = start_ptr_y_surf2 - slope_surf2_prep * start_ptr_x_surf2

    #segement (prep + end_ptr)
    b_end_prep = end_ptr_y_surf2 - slope_surf2_prep * end_ptr_x_surf2

    return slope_surf2,slope_surf2_prep,b,b_start_prep, b_end_prep

#-----------------------#
# Physical Line "surf2" #
#-----------------------#

#End Ptrs Coordinates
ptr_i_x_surf2 = 0
ptr_i_y_surf2 = 0
ptr_j_x_surf2 = 10392.3
ptr_j_y_surf2 = -6000

#-----------------------#
# Physical Line "surf4" #
#-----------------------#

#End Ptrs Coordinates
start_ptr_surf4 = -16000.0
end_ptr_surf4 = 12000.0

#-----------------------#
# Physical Line "surf6" #
#-----------------------#

#End Ptrs Coordinates
# start_ptr_surf6 = 0.0
# end_ptr_surf6 = 50000.0

#-----------------------#
# Physical Line "surf8" #
#-----------------------#

#End Ptrs Coordinates
# ptr_i_x_surf8 = -33921.924
# ptr_i_y_surf8 = -21196.771
# ptr_j_x_surf8 = 0.0
# ptr_j_y_surf8 = 0.0

#--surf2--#

#fault line slope_surf2 
slope_surf2 = (ptr_j_y_surf2 - ptr_i_y_surf2) / (ptr_j_x_surf2 - ptr_i_x_surf2)

# determine start/end ptr based on slope_surf2
if slope_surf2 < 0:
    if ( ptr_i_x_surf2 < ptr_j_x_surf2 ): #min x ptr is start_ptr
        start_ptr_x_surf2 = ptr_i_x_surf2
        start_ptr_y_surf2 = ptr_i_y_surf2
        end_ptr_x_surf2 = ptr_j_x_surf2
        end_ptr_y_surf2 = ptr_j_y_surf2
    else:
        start_ptr_x_surf2 = ptr_j_x_surf2
        start_ptr_y_surf2 = ptr_j_y_surf2
        end_ptr_x_surf2 = ptr_i_x_surf2
        end_ptr_y_surf2 = ptr_i_y_surf2
else:
    if ( ptr_i_x_surf2 < ptr_j_x_surf2 ): #max x ptr is start_ptr
        start_ptr_x_surf2 = ptr_j_x_surf2
        start_ptr_y_surf2 = ptr_j_y_surf2
        end_ptr_x_surf2 = ptr_i_x_surf2
        end_ptr_y_surf2 = ptr_i_y_surf2
    else:
        start_ptr_x_surf2 = ptr_i_x_surf2
        start_ptr_y_surf2 = ptr_i_y_surf2
        end_ptr_x_surf2 = ptr_j_x_surf2
        end_ptr_y_surf2 = ptr_j_y_surf2

slope_surf2,slope_surf2_prep,b_surf2,b_start_surf2,b_end_surf2 = find_prep_line_eq(start_ptr_x_surf2,start_ptr_y_surf2,end_ptr_x_surf2,end_ptr_y_surf2)

#---------#

# Get Physical TRIA Element Indices
# e.g. : {'surf2': {'triangle': array([    0,     1,     2, ..., 23215, 23216, 23217], dtype=uint64)}
phy_nodal_ind_surf2 = m.cell_sets_dict['surf2']['triangle']
phy_nodal_ind_surf4 = m.cell_sets_dict['surf4']['triangle']
# print(phy_nodal_ind_surf2)

# Get Crack Nodal Connectivity
phy_nodal_arr_surf2 = np.unique(np.ravel(m.cells_dict['triangle'][phy_nodal_ind_surf2,:]))
phy_nodal_arr_surf4 = np.unique(np.ravel(m.cells_dict['triangle'][phy_nodal_ind_surf4,:]))

# Find Elements Along the Crack Line
# surf2
ind_elem_crack_surf2_upper = []
ind_elem_crack_surf2_lower = []
# surf4
ind_elem_crack_surf4_upper = []
ind_elem_crack_surf4_lower = []

num_elem = np.shape(tetra_elem_connect)[0] # elem num
for elem_ind in range(num_elem):  

    print("Progress: ",round(elem_ind/num_elem*100,3),"%") 
    elem_connect_i = tetra_elem_connect[elem_ind,:] # get connectivity for current element
    
    #surf2
    num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr_surf2,elem_connect_i)) # count
    if (num_ptr_collapse == 3): # this is a crack aligned element
        coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
        coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
        #if ( len(np.argwhere(slope_surf2_prep * coord_data_x + b_start_surf2 - coord_data_y >= 0 )) > 1 and len(np.argwhere(slope_surf2_prep * coord_data_x + b_end_surf2 - coord_data_y <= 0)) > 1 ): #slope_surf2_prep > 0
        if ( slope_surf2 * np.sum(coord_data_x)/4 + b_surf2 - np.sum(coord_data_y)/4 < 0 ): #hardcode #above the fault #center ptr
            ind_elem_crack_surf2_upper.append(elem_ind)
        else:
            ind_elem_crack_surf2_lower.append(elem_ind)
        #else: 
        #    continue
    
    #surf4
    num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr_surf4,elem_connect_i)) # count
    if (num_ptr_collapse == 3): # this is a crack aligned element
        coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
        coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
        # if ( len(np.argwhere(coord_data_x > start_ptr_surf4)) >= 1 and len(np.argwhere(coord_data_x < end_ptr_surf4)) >= 1 ): #hardcode #within fault region
        if ( np.sum(coord_data_y) > 0 ): #hardcode #above the fault
            ind_elem_crack_surf4_upper.append(elem_ind)
        else:
            ind_elem_crack_surf4_lower.append(elem_ind)
        # else:
            # continue

ind_elem_all_surf2 = ind_elem_crack_surf2_upper + ind_elem_crack_surf2_lower  
ind_elem_all_surf4 = ind_elem_crack_surf4_upper + ind_elem_crack_surf4_lower

ind_subdomainid_all_surf2 = []
ind_subdomainid_all_surf4 = []
#
for i in range(len(ind_elem_crack_surf2_upper)):
    ind_subdomainid_all_surf2.append(100)
for j in range(len(ind_elem_crack_surf2_lower)):
    ind_subdomainid_all_surf2.append(200)
#
for k in range(len(ind_elem_crack_surf4_upper)):
    ind_subdomainid_all_surf4.append(300)
for l in range(len(ind_elem_crack_surf4_lower)):
    ind_subdomainid_all_surf4.append(400)

#save elem ids
np.savetxt('./list_elemid_surf2.txt',
            ind_elem_all_surf2,
            fmt='%i',
            newline=" ")
#save subdomain id
np.savetxt('./list_subdomainid_surf2.txt',
            ind_subdomainid_all_surf2,
            fmt='%i',
            newline=" ")
#
#save elem ids
np.savetxt('./list_elemid_surf4.txt',
            ind_elem_all_surf4,
            fmt='%i',
            newline=" ")
#save subdomain id
np.savetxt('./list_subdomainid_surf4.txt',
            ind_subdomainid_all_surf4,
            fmt='%i',
            newline=" ")

# print(np.shape(m.cells_dict['triangle']))
# print("--------")
# print(m.points[1610,:],"*",m.points[4,:],"*",m.points[15098,:])