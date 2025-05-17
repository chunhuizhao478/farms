#For TRIA3
import matplotlib.pyplot as plt
import numpy as np
import math
# from labellines import labelLines
import pandas as pd
import meshio
import csv

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

def stress_transformation(sigmaxx, sigmaxy, sigmayy, angle):

    #local normal stress
    local_normal_sts = 0.5 * ( sigmaxx + sigmayy ) - 0.5 * ( sigmaxx - sigmayy ) * np.cos(2*np.radians(angle)) + sigmaxy * np.sin(2*np.radians(angle))

    #local shear stress
    local_shear_sts = -0.5 * ( sigmaxx - sigmayy ) * np.sin(2*np.radians(angle)) + sigmaxy * np.cos(2*np.radians(angle))

    # print("============================================")
    # print("angle to horizontal axis: ",angle)
    # print("normal stress: ",local_normal_sts)
    # print("shear stress: ",local_shear_sts)
    # print("shear to normal ratio: ",local_shear_sts/local_normal_sts)
    # print("============================================")
    
    return local_normal_sts, local_shear_sts

def main_func(sigmayy, nu, added_normalsts, sigmayy_label, added_normalsts_label, mesh_path, branch_end_point):

    m = meshio.read(mesh_path)

    #---------------------------#
    #    Material Property      #
    #---------------------------#
    sigmaxx = nu * sigmayy #Pa #@
    sigmaxy = 0       #Pa
    sigmayy = sigmayy #Pa
    mu_s = 0.7
    mu_d = 0.1
    Dc   = 3.2e-5     #m

    #---------------------------#
    # Physical Line "embeded1"  #
    #---------------------------#

    ptr_i_x = 0
    ptr_i_y = 0.037476
    ptr_j_x = 0.141427
    ptr_j_y = 0.115869

    #fault line slope 
    slope = (ptr_j_y - ptr_i_y) / (ptr_j_x - ptr_i_x)
    radian = np.arctan(slope)
    degree = np.degrees(radian)

    #local stress field
    sigma_N, tau_S = stress_transformation(sigmaxx, sigmaxy, sigmayy, degree)

    print("============================================")
    print("angle to horizontal axis: ",degree)
    print("normal stress: ",sigma_N)
    print("shear stress: ",tau_S)
    print("shear to normal ratio: ",tau_S/sigma_N)
    print("============================================")

    #overstress 1%
    tau_S_nucleation = 1.05 * mu_s * sigma_N

    #print frictional length scale
    G = 2.17e9 #Pa

    L = G * Dc / ( ( mu_s - mu_d ) * abs(sigma_N) )
    print("frictional length scale (mm): ",L * 1e3)
    # exit(0)

    #group arr_data
    arr_data = [tau_S, sigma_N, mu_s, mu_d, Dc]

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
    else:
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

    slope,slope_prep,b_30,b_start_30,b_end_30 = find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y)

    # Get Physical Line Element Indices
    phy_nodal_ind_30 = m.cell_sets_dict['embeded1']['line']
    #print(phy_nodal_ind)

    # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
    phy_nodal_arr_30 = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind_30,:]))
    #print(phy_nodal_arr)

    # Get QUAD4 Element Connectivity
    tria_elem_connect = m.cells_dict['triangle']
    num_elem = np.shape(tria_elem_connect)[0]
    arr_elem_properties = np.zeros((num_elem,len(arr_data)))

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
                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                    arr_elem_properties[elem_ind,:] = arr_data
                    #apply nucleation shear stress
                    # x_max_lim = (61.223+10.0) * 1e-3
                    # x_min_lim = (61.223-10.0) * 1e-3
                    # y_max_lim = (71.413+5.5) * 1e-3
                    # y_min_lim = (71.413-5.5) * 1e-3
                    x_max_lim = (61.223+8.0) * 1e-3
                    x_min_lim = (61.223-8.0) * 1e-3
                    y_max_lim = (71.413+4.5) * 1e-3
                    y_min_lim = (71.413-4.5) * 1e-3
                    if ( np.sum(coord_data_x)/3 <= x_max_lim and np.sum(coord_data_x)/3 >= x_min_lim and np.sum(coord_data_y)/3 <= y_max_lim and np.sum(coord_data_y)/3 >= y_min_lim ):
                        arr_elem_properties[elem_ind,0] = tau_S_nucleation #Pa
            else:
                continue

    #only store branch
    ind_elem_all_embeded1 = ind_elem_crack_30_upper + ind_elem_crack_30_lower

    ind_subdomainid_all_embeded1 = []
    for k in range(len(ind_elem_crack_30_upper)):
        ind_subdomainid_all_embeded1.append(100)
    for l in range(len(ind_elem_crack_30_lower)):
        ind_subdomainid_all_embeded1.append(200)

    #---------------------------#
    # Physical Line "embeded2"  #
    #---------------------------#

    ptr_i_x = 0.142127
    ptr_i_y = 0.116257
    ptr_j_x = 0.203
    ptr_j_y = 0.150

    #fault line slope 
    slope = (ptr_j_y - ptr_i_y) / (ptr_j_x - ptr_i_x)
    radian = np.arctan(slope)
    degree = np.degrees(radian)

    #local stress field
    sigma_N, tau_S = stress_transformation(sigmaxx, sigmaxy, sigmayy, degree)

    #increase the normal stress by 6MPa
    sigma_N = sigma_N + added_normalsts #@

    print("============================================")
    print("angle to horizontal axis: ",degree)
    print("normal stress: ",sigma_N)
    print("shear stress: ",tau_S)
    print("shear to normal ratio: ",tau_S/sigma_N)
    print("============================================")

    #group arr_data
    arr_data = [tau_S, sigma_N, mu_s, mu_d, Dc]

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
    else:
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

    slope,slope_prep,b_30,b_start_30,b_end_30 = find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y)

    # Get Physical Line Element Indices
    phy_nodal_ind_30 = m.cell_sets_dict['embeded2']['line']
    #print(phy_nodal_ind)

    # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
    phy_nodal_arr_30 = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind_30,:]))
    #print(phy_nodal_arr)

    # Get QUAD4 Element Connectivity
    tria_elem_connect = m.cells_dict['triangle']
    # num_elem = np.shape(tria_elem_connect)[0]
    # arr_elem_properties = np.zeros((num_elem,np.shape(arr_data)[1]))

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
                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                    arr_elem_properties[elem_ind,:] = arr_data
            else:
                continue

    #only store branch
    ind_elem_all_embeded2 = ind_elem_crack_30_upper + ind_elem_crack_30_lower

    ind_subdomainid_all_embeded2 = []
    for k in range(len(ind_elem_crack_30_upper)):
        ind_subdomainid_all_embeded2.append(300)
    for l in range(len(ind_elem_crack_30_lower)):
        ind_subdomainid_all_embeded2.append(200)

    #---------------------------#
    # Physical Line "embeded3"  #
    #---------------------------#

    ptr_i_x = 0.141427
    ptr_i_y = 0.116669

    #take new ptr
    #----------------------------#
    ptr_j_x = branch_end_point[0]
    ptr_j_y = branch_end_point[1]
    #----------------------------#

    #fault line slope 
    slope = (ptr_j_y - ptr_i_y) / (ptr_j_x - ptr_i_x)
    radian = np.arctan(slope)
    degree = np.degrees(radian)

    #local stress field
    sigma_N, tau_S = stress_transformation(sigmaxx, sigmaxy, sigmayy, degree)

    print("============================================")
    print("angle to horizontal axis: ",degree)
    print("normal stress: ",sigma_N)
    print("shear stress: ",tau_S)
    print("shear to normal ratio: ",tau_S/sigma_N)
    print("============================================")

    #group arr_data
    arr_data = [tau_S, sigma_N, mu_s, mu_d, Dc]

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
    else:
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

    slope,slope_prep,b_30,b_start_30,b_end_30 = find_prep_line_eq(start_ptr_x,start_ptr_y,end_ptr_x,end_ptr_y)

    # Get Physical Line Element Indices
    phy_nodal_ind_30 = m.cell_sets_dict['embeded3']['line']
    #print(phy_nodal_ind)

    # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
    phy_nodal_arr_30 = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind_30,:]))
    #print(phy_nodal_arr)

    # Get QUAD4 Element Connectivity
    tria_elem_connect = m.cells_dict['triangle']
    # num_elem = np.shape(tria_elem_connect)[0]
    # arr_elem_properties = np.zeros((num_elem,np.shape(arr_data)[1]))

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
                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                    arr_elem_properties[elem_ind,:] = arr_data
            else:
                continue

    #only store branch
    ind_elem_all_embeded3 = ind_elem_crack_30_upper + ind_elem_crack_30_lower

    ind_subdomainid_all_embeded3 = []
    for k in range(len(ind_elem_crack_30_upper)):
        ind_subdomainid_all_embeded3.append(300)
    for l in range(len(ind_elem_crack_30_lower)):
        ind_subdomainid_all_embeded3.append(100)

    ind_elem_all_embeded = ind_elem_all_embeded1 + ind_elem_all_embeded2 + ind_elem_all_embeded3
    ind_subdomainid_all_embeded = ind_subdomainid_all_embeded1 + ind_subdomainid_all_embeded2 + ind_subdomainid_all_embeded3

    #save elem ids
    np.savetxt("./case_P"+str(sigmayy_label)+"/case_deltasigma"+str(added_normalsts_label)+"/list_elemid_embeded.txt",
                ind_elem_all_embeded,
                fmt='%i',
                newline=" ")

    #save subdomain id
    np.savetxt("./case_P"+str(sigmayy_label)+"/case_deltasigma"+str(added_normalsts_label)+"/list_subdomainid_embeded.txt",
                ind_subdomainid_all_embeded,
                fmt='%i',
                newline=" ")

    #save property elem csv
    with open("./case_P"+str(sigmayy_label)+"/case_deltasigma"+str(added_normalsts_label)+"/faults_elem_properties.csv", 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)

        # write multiple rows
        writer.writerows(arr_elem_properties)

# exit(0)

