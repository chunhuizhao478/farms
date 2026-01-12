#For TRIA3

import numpy as np
import meshio
import os
import shutil

m = meshio.read('tpv142d_100m.msh')

subfolder = "raw_results"

if os.path.exists(subfolder):
    # Remove all files and subfolders inside the "result" folder.
    for filename in os.listdir(subfolder):
        file_path = os.path.join(subfolder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)  # Remove a file or a symbolic link.
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)  # Remove a directory and all its contents.
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")
else:
    # Create the subfolder if it does not exist.
    os.makedirs(subfolder)

#m.cells_dict['line']

# print(m.cells_dict)

#---------------------#
# Physical Line "180" #
#---------------------#

# Start/End
start_ptr = -16000
end_ptr = 12000

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
# Save the coordinates of the crack line
xcoord_crack_10 = []
ycoord_crack_10 = []
num_elem = np.shape(tria_elem_connect)[0] # elem num
for elem_ind in range(num_elem):   
    elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
    #print(elem_connect_i)
    num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr,elem_connect_i)) # count
    if (num_ptr_collapse >= 1): # this is a crack aligned element
        coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
        coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
        if ( len(np.argwhere(coord_data_x > start_ptr)) >= 1 and len(np.argwhere(coord_data_x < end_ptr)) >= 1 ): #hardcode #within fault region
            if ( np.sum(coord_data_y) > 0 ): #hardcode #above the fault
                ind_elem_crack_10_upper.append(elem_ind)
                # Save the coordinates of the crack line
                xcoord_crack_10.append(coord_data_x)
                ycoord_crack_10.append(coord_data_y)
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

np.savetxt(subfolder+'/res_xcoord_crack_mf.txt', xcoord_crack_10, fmt='%s')
np.savetxt(subfolder+'/res_ycoord_crack_mf.txt', ycoord_crack_10, fmt='%s')
np.savetxt(subfolder+'/res_ind_elem_crack_mf.txt', ind_elem_crack_10_upper, fmt='%s')

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
start_ptr_x = 86.6025
start_ptr_y = -50

end_ptr_x = 10392.3
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
# Save the coordinates of the crack line
xcoord_crack_30 = []
ycoord_crack_30 = []
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
                # Save the coordinates of the crack line
                xcoord_crack_30.append(coord_data_x)
                ycoord_crack_30.append(coord_data_y)
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

np.savetxt(subfolder+'/res_xcoord_crack_bf.txt', xcoord_crack_30, fmt='%s')
np.savetxt(subfolder+'/res_ycoord_crack_bf.txt', ycoord_crack_30, fmt='%s')
np.savetxt(subfolder+'/res_ind_elem_crack_bf.txt', ind_elem_crack_30_upper, fmt='%s')

exit(0)