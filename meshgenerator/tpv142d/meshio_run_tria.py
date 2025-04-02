#For TRIA3

import numpy as np
import meshio

m = meshio.read('tpv142d_100m.msh')

#m.cells_dict['line']

# print(m.cells_dict)

#---------------------#
# Physical Line "180" #
#---------------------#

# Start/End
start_ptr = -10000
end_ptr = 18000

# Get Physical Line Element Indices
phy_nodal_ind = m.cell_sets_dict['180']['line']
#print(phy_nodal_ind)

# Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
phy_nodal_arr = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind,:]))
#print(phy_nodal_arr)

# Get QUAD4 Element Connectivity #Need to define Physical Surface
tria_elem_connect = m.cells_dict['triangle']

# Find Elements Along the Crack Line
ind_elem_crack_10_upper = []
ind_elem_crack_10_lower = []
num_elem = np.shape(tria_elem_connect)[0] # elem num
for elem_ind in range(num_elem):   
    elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
    #print(elem_connect_i)
    num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr,elem_connect_i)) # count
    if (num_ptr_collapse >= 1): # this is a crack aligned element
        coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
        coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
        if ( len(np.argwhere(coord_data_x >= start_ptr)) >= 1 and len(np.argwhere(coord_data_x <= end_ptr)) >= 1 ): #hardcode #within fault region
            if ( np.sum(coord_data_y) > 0 ): #hardcode #above the fault
                ind_elem_crack_10_upper.append(elem_ind)
            else:
                ind_elem_crack_10_lower.append(elem_ind)
        else:
            continue

#output
print("ind_elem_crack_180")
for item in ind_elem_crack_10_upper:
    print(item,end=" ")
for item2 in ind_elem_crack_10_lower:
    print(item2,end=" ")

print("\n ind_elem_crack_180[subdomainID]")
for i in range(len(ind_elem_crack_10_upper)):
    print(100,end=" ")
for j in range(len(ind_elem_crack_10_lower)):
    print(200,end=" ")
print("\n 180 node number: ",len(ind_elem_crack_10_lower))


def find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y):

    #fault line slope 
    slope = (end_ptr_y - start_ptr_y) / (end_ptr_x - start_ptr_x)

    #prep line slope
    slope_prep = - 1 / slope

    #segement 
    b = start_ptr_y - slope * start_ptr_x

    #segement (prep + start_ptr)
    b_start_prep = start_ptr_y - slope_prep * start_ptr_x

    #segement (prep + end_ptr)
    b_end_prep = end_ptr_y - slope_prep * end_ptr_x

    return slope,slope_prep,b,b_start_prep, b_end_prep

#---------------------#
# Physical Line "30"  #
#---------------------#

# Start/End
start_ptr_x = 6086.6025
start_ptr_y = -50

end_ptr_x = 16392.3
end_ptr_y = -6000

slope,slope_prep,b_30,b_start_30,b_end_30 = find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y)

# Get Physical Line Element Indices
phy_nodal_ind_30 = m.cell_sets_dict['30']['line']
#print(phy_nodal_ind)

# Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
phy_nodal_arr_30 = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind_30,:]))
#print(phy_nodal_arr)

# Get QUAD4 Element Connectivity
tria_elem_connect = m.cells_dict['triangle']

# Find Elements Along the Crack Line
ind_elem_crack_30_upper = []
ind_elem_crack_30_lower = []
num_elem = np.shape(tria_elem_connect)[0] # elem num
for elem_ind in range(num_elem):   
    elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
    num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr_30,elem_connect_i)) # count
    if (num_ptr_collapse >= 1): # this is a crack aligned element
        coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
        coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
        if ( len(np.argwhere(slope_prep * coord_data_x + b_start_30 - coord_data_y >= 0 )) > 1 and len(np.argwhere(slope_prep * coord_data_x + b_end_30 - coord_data_y <= 0)) > 1 ): #slope_prep > 0
            if ( slope * np.sum(coord_data_x)/3 + b_30 - np.sum(coord_data_y)/3 < 0 ): #hardcode #above the fault #center ptr
                ind_elem_crack_30_upper.append(elem_ind)
            else:
                ind_elem_crack_30_lower.append(elem_ind)
        else:
            continue

#output
print("\n ind_elem_crack_30")
for item in ind_elem_crack_30_upper:
    print(item,end=" ")
for item2 in ind_elem_crack_30_lower:
    print(item2,end=" ")

print("\n ind_elem_crack_30[subdomainID]")
for i in range(len(ind_elem_crack_30_upper)):
    print(300,end=" ")
for j in range(len(ind_elem_crack_30_lower)):
    print(400,end=" ")
print("\n 30 node number: ",len(ind_elem_crack_30_lower))

exit(0)