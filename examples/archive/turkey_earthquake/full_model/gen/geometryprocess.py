#scale of maps 1:1251654
#read vertices locations
import numpy as np
import math
import matplotlib.pyplot as plt
from labellines import labelLines
import csv
import gmsh
import meshio
import pandas as pd

def helper_rotate(coordx,coordy,angle):

    # Define angle in radians
    theta = math.radians(angle)

    # Define rotation matrix
    R = [[math.cos(theta), -math.sin(theta)],
        [math.sin(theta), math.cos(theta)]]

    # Multiply original coordinates by rotation matrix
    x_new = np.multiply(R[0][0],coordx) + np.multiply(R[0][1],coordy)
    y_new = np.multiply(R[1][0],coordx) + np.multiply(R[1][1],coordy)

    # plt.figure()
    # plt.plot(coordx,coordy,"b*")
    # # plt.plot(x_new,y_new,"r*")
    # plt.title("Fault Trace Geometry")
    # plt.xlabel("x coordinate")
    # plt.xlabel("y coordinate")
    # plt.show()
    
    # plt.figure()
    # # plt.plot(coordx,coordy,"b*")
    # plt.plot(x_new,y_new,"r*")
    # plt.title("Fault Trace Geometry")
    # plt.xlabel("x coordinate")
    # plt.xlabel("y coordinate")
    # plt.show()

    return x_new, y_new

def helper_physical(coordx,coordy,shriftdist_x,shriftdist_y,scale):

    #shrift origin
    coordxnew = coordx - shriftdist_x
    coordynew = coordy - shriftdist_y

    #physical coordinate
    coordxnew = np.multiply(coordxnew,scale)
    coordynew = np.multiply(coordynew,scale)

    # plt.figure()
    # # plt.plot(coordx,coordy,"b*")
    # plt.plot(coordxnew,coordynew,"r*")
    # plt.title("Fault Trace Geometry")
    # plt.xlabel("x coordinate")
    # plt.xlabel("y coordinate")
    # plt.show()

    return coordxnew,coordynew

def helper_movecoords(coordx,coordy):


    coordx_new = np.zeros(np.shape(coordx))
    coordy_new = np.zeros(np.shape(coordx))

    #make junction part main fault striaght
    moveindex = []
    for i in range(len(coordx)):
        if coordx[i] < -2.5e4 and coordy[i] < 10000:
            coordx_new[i] = coordx[i] + abs(3e4-5.2057e4)
            coordy_new[i] = coordy[i] - abs(4.2e3) 
            moveindex.append(i)
        else:
            coordx_new[i] = coordx[i]
            coordy_new[i] = coordy[i]

    # plt.figure()
    # plt.plot(coordx_new[moveindex],coordy_new[moveindex],"b.")
    # plt.plot(coordx,coordy,"r.")
    # plt.title("Fault Trace Geometry (move step 1)")
    # plt.xlabel("x coordinate")
    # plt.ylabel("y coordinate")
    # plt.show()
    
    #move south end
    for i in range(len(coordx_new)):
        if coordx_new[i] < -8.38e4 and coordy_new[i] < -3.75e4:
            coordx_new[i] = coordx_new[i] - abs(17000)
            coordy_new[i] = coordy_new[i] - abs(12261) 
            moveindex.append(i)
        else:
            coordx_new[i] = coordx_new[i]
            coordy_new[i] = coordy_new[i]

    # plt.figure()
    # plt.plot(coordx_new[moveindex],coordy_new[moveindex],"b.")
    # plt.plot(coordx,coordy,"r.")
    # plt.title("Fault Trace Geometry (move step 2)")
    # plt.xlabel("x coordinate")
    # plt.ylabel("y coordinate")
    # plt.show()

    return coordx_new,coordy_new

def helper_addlines(coordxm,coordym,dictlines):

    plt.figure()
    plt.plot(coordxm,coordym,"r*")
    for i in range(1,len(dictlines)+1):
        print([dictlines[i][0],dictlines[i][2]])
        print([dictlines[i][1],dictlines[i][3]])
        plt.plot([dictlines[i][0],dictlines[i][2]],[dictlines[i][1],dictlines[i][3]],"b-",label=str(i))
    labelLines(plt.gca().get_lines(), zorder=2.5,yoffsets=1)
    plt.title("Fault Trace Geometry")
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.show()

def helper_getallptrs(dictlines,mesh_type):

    list_locs = []
    for i in range(1,len(dictlines)+1):
        list_locs.extend(dictlines[i])
    
    arr_locs = np.array(list_locs).reshape(int(len(list_locs)/2),2)

    arr_unique_locs = np.unique(arr_locs, axis=0)

    dict_unique_ptrs = {}

    # plt.figure()
    for i in range(np.shape(arr_unique_locs)[0]):
        # plt.plot(arr_unique_locs[i][0],arr_unique_locs[i][1],'r*')
        # plt.annotate(str(i+1),  # this is the text (put lab here to use tlab as string)
                #  (arr_unique_locs[i][0], arr_unique_locs[i][1]),  # this is the point to label
                #  textcoords="offset points",  # how to position the text
                #  xytext=(0, 5),  # distance from text to points (x,y)
                #  ha='center')
        dict_unique_ptrs[i+1] = [arr_unique_locs[i][0],arr_unique_locs[i][1]]
    # plt.title("Fault Piecewise Constant Geometry Kink Points")
    # plt.xlabel("x coordinate")
    # plt.ylabel("y coordinate")
    # plt.show()

    return dict_unique_ptrs

def helper_addmeshlines(dict_lines_connectivity,dict_unique_ptrs,coordmx,coordmy):

    # plt.figure()
    for i in dict_unique_ptrs.keys():
        plt.plot(dict_unique_ptrs[i][0],dict_unique_ptrs[i][1],'b.')
        plt.annotate(str(i),  # this is the text (put lab here to use tlab as string)
                 (dict_unique_ptrs[i][0], dict_unique_ptrs[i][1]),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 5),  # distance from text to points (x,y)
                 ha='center')
    # plt.plot(coordxm,coordym,"r*")
    # plt.show()
    # for i in range(1,13+1):
    for i in dict_lines_connectivity.keys():
        line_i = dict_lines_connectivity[i]
        line_i_start = dict_unique_ptrs[line_i[0]]
        line_i_end = dict_unique_ptrs[line_i[1]]
        plt.plot([line_i_start[0],line_i_end[0]],[line_i_start[1],line_i_end[1]],"b-",label=str(i))
    labelLines(plt.gca().get_lines(), zorder=2.5,yoffsets=1)
    plt.title("Fault Trace Geometry (addmeshlines check)")
    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    plt.show()

def helper_reverse(coordxp,coordyp,dict_unique_ptrs,dict_lines_connectivity,scale,shriftdist_x,shriftdist_y,angle):

    #help reverse the process to find model coordinates on the map

    #divide by the scale
    coordxp_scaleback = np.multiply(coordxp,1/scale)
    coordyp_scaleback = np.multiply(coordyp,1/scale)
    for i in dict_unique_ptrs.keys():
        dict_unique_ptrs[i][0] = dict_unique_ptrs[i][0] * 1/scale
        dict_unique_ptrs[i][1] = dict_unique_ptrs[i][1] * 1/scale

    #shrift back
    coordxp_shriftback = coordxp_scaleback + shriftdist_x
    coordyp_shriftback = coordyp_scaleback + shriftdist_y
    for i in dict_unique_ptrs.keys():
        dict_unique_ptrs[i][0] = dict_unique_ptrs[i][0] + shriftdist_x
        dict_unique_ptrs[i][1] = dict_unique_ptrs[i][1] + shriftdist_y

    #arr_unique_ptrs
    arr_unique_ptrs_x = np.zeros((len(dict_unique_ptrs),1))
    arr_unique_ptrs_y = np.zeros((len(dict_unique_ptrs),1))

    count = 0
    for i in dict_unique_ptrs.keys():
        arr_unique_ptrs_x[count] = dict_unique_ptrs[i][0] 
        arr_unique_ptrs_y[count] = dict_unique_ptrs[i][1] 
        count = count + 1

    # plt.figure(figsize=(20,10))
    # plt.plot(coordxp_shriftback,coordyp_shriftback,"r*")
    # plt.plot(arr_unique_ptrs_x,arr_unique_ptrs_y,"b.")
    # plt.show()

    #rotate back
    angle = -1 * angle

    coordxp_rotateback,coordyp_rotateback = helper_rotate(coordxp_shriftback,coordyp_shriftback,angle)
    uniqueptrsx_rotateback,uniqueptrsy_rotateback = helper_rotate(arr_unique_ptrs_x,arr_unique_ptrs_y,angle)

    #put back to dict
    count = 0
    for i in dict_unique_ptrs.keys():
        dict_unique_ptrs[i][0] = uniqueptrsx_rotateback[count]
        dict_unique_ptrs[i][1] = uniqueptrsy_rotateback[count] 
        count = count + 1

    #adjustment
    for i in dict_unique_ptrs.keys():
        dict_unique_ptrs[i][0] = dict_unique_ptrs[i][0]
        dict_unique_ptrs[i][1] = dict_unique_ptrs[i][1] 

    plt.figure(figsize=(20,10))
    for i in dict_unique_ptrs.keys():
        plt.plot(dict_unique_ptrs[i][0],dict_unique_ptrs[i][1],'b.')
        plt.annotate(str(i),  # this is the text (put lab here to use tlab as string)
                 (dict_unique_ptrs[i][0], dict_unique_ptrs[i][1]),  # this is the point to label
                 textcoords="offset points",  # how to position the text
                 xytext=(0, 5),  # distance from text to points (x,y)
                 ha='center')
    for i in dict_lines_connectivity.keys():
        line_i = dict_lines_connectivity[i]
        line_i_start = dict_unique_ptrs[line_i[0]]
        line_i_end = dict_unique_ptrs[line_i[1]]
        plt.plot([line_i_start[0],line_i_end[0]],[line_i_start[1],line_i_end[1]],"b-",label=str(i))
    plt.plot(coordxp_rotateback,coordyp_rotateback,"r*")
    plt.show()

    #save points
    iden_col = np.linspace(1,len(dict_unique_ptrs),len(dict_unique_ptrs)).reshape(np.shape(uniqueptrsx_rotateback))
    arr_cols = np.hstack((uniqueptrsx_rotateback,uniqueptrsy_rotateback))
    arr_cols = np.hstack((arr_cols,iden_col))
    header = ["Lat","Long","Index"]
    with open('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/model_vertices_map.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # write the header
            writer.writerow(header)

            # write multiple rows
            writer.writerows(arr_cols)



def helper_Hangle(dict_lines_connectivity,dict_unique_ptrs):

    #calcuate each line angle in degree with respect to horizontal line
    #return a list of angles

    #initialize numpy array for storing 
    arr_lines_angles = np.zeros((len(dict_lines_connectivity),1))

    #loop over lines
    for i in range(1,len(dict_lines_connectivity)+1):
        line_i = dict_lines_connectivity[i]
        line_i_start = dict_unique_ptrs[line_i[0]]
        line_i_end = dict_unique_ptrs[line_i[1]]
        slope = (line_i_end[1]-line_i_start[1])/(line_i_end[0]-line_i_start[0])
        radian = np.arctan(slope)
        degree = np.degrees(radian)
        arr_lines_angles[i-1] = degree

    print(arr_lines_angles)

    return arr_lines_angles

def propertycalc_sts(psts1,psts3,arr_lines_angles,angleofmaxcomp):

    #inherit from propcalcv1.py in propertycalc folder

    #initialize containers for normal stress, shear stress, shear-normal ratio
    arr_normalsts = np.zeros(np.shape(arr_lines_angles))
    arr_shearsts = np.zeros(np.shape(arr_lines_angles))
    arr_sheartonormal = np.zeros(np.shape(arr_lines_angles))
    
    #maximum compression angle wrt each fault
    #take absolute value
    arr_angleofmaxcompwfault = angleofmaxcomp - arr_lines_angles

    ##Note:
    #clockwise to be positive
    #counterclockwise to be negative
    
    #loop over lines
    for i in range(np.shape(arr_angleofmaxcompwfault)[0]):

        #current angle
        theta = arr_angleofmaxcompwfault[i]

        #define cost sint
        cost = math.cos(math.radians(theta))
        sint = math.sin(math.radians(theta))

        #normal stress
        sigma_n = psts1 * sint ** 2 + psts3 * cost ** 2
        arr_normalsts[i] = sigma_n

        #shear stress
        tau_s = ( psts1 - psts3 ) * cost * sint
        arr_shearsts[i] = tau_s

        #check ratio #this is absolute quantity
        ratio_shear_normal = ( 1 - psts3 / psts1 ) * cost * sint / ( sint ** 2 + ( psts3 / psts1 ) * cost ** 2 )
        arr_sheartonormal[i] = abs(ratio_shear_normal)

        print("============================================")
        print("segement: ",i+1)
        print("angle to horizontal axis: ",arr_lines_angles[i])
        print("angle to max compression: ",arr_angleofmaxcompwfault[i])
        print("normal stress: ",arr_normalsts[i])
        print("shear stress: ",arr_shearsts[i])
        print("shear to normal ratio: ",arr_sheartonormal[i])
        print("============================================")

    # print("mf check: ")
    # print("normal stress: ",arr_normalsts[30-1])
    # print("shear stress: ",arr_shearsts[30-1])
    # print("splay check: ")
    # print("normal stress: ",arr_normalsts[31-1])
    # print("shear stress: ",arr_shearsts[31-1])

    return arr_angleofmaxcompwfault, arr_normalsts, arr_shearsts, arr_sheartonormal

def propertycalc_friction(arr_normalsts, 
                          arr_Sparam, 
                          arr_sheartonormal, 
                          ind_importantsegments_mainfaultbykink,
                          ind_mainfault,
                          ind_splayfault,
                          ind_mainfaultbeyondkink):

    #compute friction property
    #tentative mu_s
    arr_tetative_mus = np.ones(np.shape(arr_angleofmaxcompwfault))

    #constraint 1: mu_s > mu on main faults beyond kink
    max_mu_beyondkink = max(arr_sheartonormal[ind_importantsegments_mainfaultbykink])
    max_mu_lineseg = np.argmax(arr_sheartonormal[ind_importantsegments_mainfaultbykink])

    print("maximum mu beyond kink value: ", max_mu_beyondkink)
    print("maximum mu beyond kink line segement: ", ind_importantsegments_mainfaultbykink[ind_mainfaultbeyondkink[max_mu_lineseg]])

    #pick a fix value + 0.1
    diff = 0.1
    max_mus_beyond_kink = round(max_mu_beyondkink[0] + diff)

    arr_tetative_mus = arr_tetative_mus * max_mus_beyond_kink
    
    print("choose mu_s beyond kink: ",max_mus_beyond_kink)

    #tentative mu_d
    arr_tetative_mud = np.zeros(np.shape(arr_angleofmaxcompwfault))

    return


def gmsh_add_points(dict_unique_ptrs,lc,mesh_type):

    gmsh.initialize()

    if mesh_type == "full":
    #add embed points -> boundary points
        for i in range(1,len(dict_unique_ptrs)+1):
            node_i = dict_unique_ptrs[i]
            node_i_xcoord = node_i[0]
            node_i_ycoord = node_i[1]
            gmsh.model.geo.addPoint(node_i_xcoord,node_i_ycoord,0,lc,tag=i)
    else:
        for i in dict_unique_ptrs.keys():
            node_i = dict_unique_ptrs[i]
            node_i_xcoord = node_i[0]
            node_i_ycoord = node_i[1]
            gmsh.model.geo.addPoint(node_i_xcoord,node_i_ycoord,0,lc,tag=i)

def gmsh_add_elems(dict_lines_connectivity,mesh_type):

    if mesh_type == "full":
        for i in range(1,len(dict_lines_connectivity)+1):
            ptr_i = dict_lines_connectivity[i][0]
            ptr_j = dict_lines_connectivity[i][1]
            gmsh.model.geo.addLine(ptr_i,ptr_j,tag=i)
    else:
        for i in dict_lines_connectivity.keys():
            ptr_i = dict_lines_connectivity[i][0]
            ptr_j = dict_lines_connectivity[i][1]
            gmsh.model.geo.addLine(ptr_i,ptr_j,tag=i)

def gmsh_add_surfaces(line_loop):

    loopi = gmsh.model.geo.addCurveLoop(line_loop,tag=1)

    gmsh.model.geo.addPlaneSurface([loopi],tag=1)

def gmsh_embed_entities(dict_unique_ptrs,dict_lines_connectivity,mesh_type):

    #We have to synchronize beforce embedding entities:
    gmsh.model.geo.synchronize()

    if mesh_type == "full":
        #embed points
        for i in range(1,len(dict_unique_ptrs)+1-4):
            gmsh.model.mesh.embed(0,[i],2,1)
        
        #embed lines
        for i in range(1,len(dict_lines_connectivity)+1-4):
            gmsh.model.mesh.embed(1,[i],2,1)
    elif mesh_type == "leftpostkink": #hardcode
        #embed points
        for i in range(1,len(dict_unique_ptrs)+1-4):
            gmsh.model.mesh.embed(0,[i],2,1)
        #embed lines
        for i in range(1,len(dict_lines_connectivity)+1-4):
            gmsh.model.mesh.embed(1,[i],2,1)
    elif mesh_type == "rightpostkink":
        #embed points
        for i in [21,22,23,24,25,26,27]: #hardcode
            gmsh.model.mesh.embed(0,[i],2,1)
        #embed lines
        for i in [21,22,23,24,25,26]: #hardcode
            gmsh.model.mesh.embed(1,[i],2,1)

    #synchronize again after embedding ptrs & elems
    gmsh.model.geo.synchronize()

def gmsh_add_physical_entities(dict_lines_connectivity, mesh_type):

    #lines
    if mesh_type == "full":
        for i in range(1,len(dict_lines_connectivity)+1-4):
            gmsh.model.addPhysicalGroup(1,[i],name="segment" + str(i))
    elif mesh_type == "leftpostkink":
        for i in range(1,len(dict_lines_connectivity)+1-4):
            gmsh.model.addPhysicalGroup(1,[i],name="segment" + str(i))
    elif mesh_type == "rightpostkink":
        for i in [21,22,23,24,25,26]:
            gmsh.model.addPhysicalGroup(1,[i],name="segment" + str(i))
    
    #surfaces
    gmsh.model.addPhysicalGroup(2,[1],name="surf" + str(100))

    #synchronize
    gmsh.model.geo.synchronize()

def gmsh_meshing(num_algo_meshing, file_path):

    gmsh.option.setNumber("Mesh.Algorithm", num_algo_meshing)

    #generate mesh
    gmsh.model.mesh.generate(dim=2)
    gmsh.write(file_path)
    gmsh.fltk.run()
    gmsh.finalize()

def meshio_subdomainids(dict_lines_connectivity,
                        dict_unique_ptrs,
                        file_path,
                        main_path_segmentids,
                        branch_path_segmentids,
                        daughter_path_segmentids,
                        arr_data,
                        splayfault_nucleation_shear_sts,
                        header):

    #read file
    m = meshio.read(file_path)

    #initialize container
    ind_elem_crack_mainpath_upper = {}
    ind_elem_crack_mainpath_lower = {}
    ind_elem_crack_branchpath_upper = {}
    ind_elem_crack_branchpath_lower = {}
    ind_elem_crack_daughterpath_upper = {}
    ind_elem_crack_daughterpath_lower = {}
    
    tria_elem_connect = m.cells_dict['triangle']
    num_elem = np.shape(tria_elem_connect)[0]
    arr_elem_properties = np.zeros((num_elem,np.shape(arr_data)[1]))

    #first main path, then branch path, then daughter path
    segmentids_sequence = main_path_segmentids + branch_path_segmentids + daughter_path_segmentids

    #arr_data SI unit
    arr_data[:,0] = arr_data[:,0] * 1e6 #Pa shear stress
    arr_data[:,1] = arr_data[:,1] * 1e6 #Pa normal stress

    #loop over lines
    for i in segmentids_sequence:

        print("current segment: ", i)

        if i in main_path_segmentids:
            ind_elem_crack_mainpath_upper[i] = []
            ind_elem_crack_mainpath_lower[i] = []
        elif i in branch_path_segmentids:
            ind_elem_crack_branchpath_upper[i] = []
            ind_elem_crack_branchpath_lower[i] = []
        else:
            ind_elem_crack_daughterpath_upper[i] = []
            ind_elem_crack_daughterpath_lower[i] = []
        
        line_i = dict_lines_connectivity[i]
        line_i_start = dict_unique_ptrs[line_i[0]]
        line_i_end = dict_unique_ptrs[line_i[1]]

        #only main fault (slope == 0)
        if i == 20 or i == 22:

            start_ptr = line_i_start[0]
            end_ptr = line_i_end[0]

            # Get Physical Line Element Indices
            phy_nodal_ind = m.cell_sets_dict['segment'+str(i)]['line']
            #print(phy_nodal_ind)

            # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
            phy_nodal_arr = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind,:]))
            #print(phy_nodal_arr)

            # Get QUAD4 Element Connectivity #Need to define Physical Surface
            tria_elem_connect = m.cells_dict['triangle']

            # Find Elements Along the Crack Line
            num_elem = np.shape(tria_elem_connect)[0] # elem num
            for elem_ind in range(num_elem):   
                elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
                #print(elem_connect_i)
                num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr,elem_connect_i)) # count
                if (num_ptr_collapse >= 1): # this is a crack aligned element
                    coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
                    coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
                    if ( len(np.argwhere(coord_data_x >= start_ptr)) > 1 and len(np.argwhere(coord_data_x < end_ptr)) > 1 ): #hardcode #within fault region
                        if ( np.sum(coord_data_y) > 0 ): #hardcode #above the fault
                            if i in main_path_segmentids:
                                ind_elem_crack_mainpath_upper[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            elif i in branch_path_segmentids:
                                ind_elem_crack_branchpath_upper[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0: #avoid overwrite
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            elif i in daughter_path_segmentids:
                                ind_elem_crack_daughterpath_upper[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            else:
                                print("ERROR! CAN NOT FIND GROUP! segement: ", i)
                                exit()
                        else:
                            if i in main_path_segmentids:
                                ind_elem_crack_mainpath_lower[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            elif i in branch_path_segmentids:
                                ind_elem_crack_branchpath_lower[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0: #avoid overwrite
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            elif i in daughter_path_segmentids:
                                ind_elem_crack_daughterpath_lower[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0: #avoid overwrite
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            else:
                                print("ERROR! CAN NOT FIND GROUP! segement: ", i)
                                exit()
                    else:
                        continue
        
        else:

            ptr_i_x = line_i_start[0]
            ptr_i_y = line_i_start[1]
            ptr_j_x = line_i_end[0]
            ptr_j_y = line_i_end[1]

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
            phy_nodal_ind_30 = m.cell_sets_dict['segment'+str(i)]['line']
            #print(phy_nodal_ind)

            # Get Crack Nodal Connectivity # plus 1 to retrieve the same index in gmsh
            phy_nodal_arr_30 = np.unique(np.ravel(m.cells_dict['line'][phy_nodal_ind_30,:]))
            #print(phy_nodal_arr)

            # Get QUAD4 Element Connectivity
            tria_elem_connect = m.cells_dict['triangle']

            # Find Elements Along the Crack Line
            num_elem = np.shape(tria_elem_connect)[0] # elem num
            for elem_ind in range(num_elem):   
                elem_connect_i = tria_elem_connect[elem_ind,:] # get connectivity for current element
                num_ptr_collapse = np.count_nonzero(np.isin(phy_nodal_arr_30,elem_connect_i)) # count
                if (num_ptr_collapse >= 1): # this is a crack aligned element
                    coord_data_x = m.points[elem_connect_i][:,0] #numpy.ndarray
                    coord_data_y = m.points[elem_connect_i][:,1] #numpy.ndarray
                    if ( len(np.argwhere(slope_prep * coord_data_x + b_start_30 - coord_data_y >= 0 )) > 1 and len(np.argwhere(slope_prep * coord_data_x + b_end_30 - coord_data_y <= 0)) > 1 ): #slope_prep > 0
                        if ( slope * np.sum(coord_data_x)/3 + b_30 - np.sum(coord_data_y)/3 < 0 ): #hardcode #above the fault #center ptr
                            if i in main_path_segmentids:
                                ind_elem_crack_mainpath_upper[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            elif i in branch_path_segmentids:
                                ind_elem_crack_branchpath_upper[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0: #avoid overwrite
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                                    #nucleation region #know splay fault in branch_path_segments
                                    if (np.sum(coord_data_x)/3 < -24763.005 and np.sum(coord_data_x)/3 > -26119.881 and np.sum(coord_data_y)/3 < -15473.643 and np.sum(coord_data_y)/3 > -16321.513):
                                        arr_elem_properties[elem_ind,0] = splayfault_nucleation_shear_sts * 1e6 #Pa
                            elif i in daughter_path_segmentids:
                                ind_elem_crack_daughterpath_upper[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            else:
                                print("ERROR! CAN NOT FIND GROUP! segement: ", i)
                                exit()
                        else:
                            if i in main_path_segmentids:
                                ind_elem_crack_mainpath_lower[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            elif i in branch_path_segmentids:
                                ind_elem_crack_branchpath_lower[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0: #avoid overwrite
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                                    #nucleation region #know splay fault in branch_path_segments
                                    if (np.sum(coord_data_x)/3 < -24763.005 and np.sum(coord_data_x)/3 > -26119.881 and np.sum(coord_data_y)/3 < -15473.643 and np.sum(coord_data_y)/3 > -16321.513):
                                        arr_elem_properties[elem_ind,0] = splayfault_nucleation_shear_sts * 1e6 #Pa
                            elif i in daughter_path_segmentids:
                                ind_elem_crack_daughterpath_lower[i].append(elem_ind)
                                if np.count_nonzero(arr_elem_properties[elem_ind,:]) == 0:
                                    arr_elem_properties[elem_ind,:] = arr_data[i-1,:]
                            else:
                                print("ERROR! CAN NOT FIND GROUP! segement: ", i)
                                exit()
                    else:
                        continue

    #sum over all
    ind_elem_mainpath = []
    ind_subdomainid_mainpath = []
    ind_elem_branchpath = []
    ind_subdomainid_branchpath = []
    ind_elem_daughterpath = []
    ind_subdomainid_daughterpath = []

    #loop
    count = 100
    for j in main_path_segmentids:
        print("create main path subdomainid for segment: ",j)
        ind_elem_mainpath.extend(ind_elem_crack_mainpath_upper[j])
        ind_subdomainid_mainpath.extend([count] * len(ind_elem_crack_mainpath_upper[j]))
        ind_elem_mainpath.extend(ind_elem_crack_mainpath_lower[j])
        ind_subdomainid_mainpath.extend([count+100] * len(ind_elem_crack_mainpath_lower[j]))
        count += 200
    for j in branch_path_segmentids:
        print("create branch path subdomainid for segment: ",j)
        ind_elem_branchpath.extend(ind_elem_crack_branchpath_upper[j])
        ind_subdomainid_branchpath.extend([count] * len(ind_elem_crack_branchpath_upper[j]))
        ind_elem_branchpath.extend(ind_elem_crack_branchpath_lower[j])
        ind_subdomainid_branchpath.extend([count+100] * len(ind_elem_crack_branchpath_lower[j]))
        count += 200
    for j in daughter_path_segmentids:
        print("create branch path subdomainid for segment: ",j)
        ind_elem_daughterpath.extend(ind_elem_crack_daughterpath_upper[j])
        ind_subdomainid_daughterpath.extend([count] * len(ind_elem_crack_daughterpath_upper[j]))
        ind_elem_daughterpath.extend(ind_elem_crack_daughterpath_lower[j])
        ind_subdomainid_daughterpath.extend([count+100] * len(ind_elem_crack_daughterpath_lower[j]))
        count += 200

    if mesh_type == "full":
        #save elem ids
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_mainpath.txt',
                    ind_elem_mainpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_branchpath.txt',
                    ind_elem_branchpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_daughterpath.txt',
                    ind_elem_daughterpath,
                    fmt='%i',
                    newline=" ")

        #save subdomain id
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_mainpath.txt',
                    ind_subdomainid_mainpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_branchpath.txt',
                    ind_subdomainid_branchpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_daughterpath.txt',
                    ind_subdomainid_daughterpath,
                    fmt='%i',
                    newline=" ")

        #save property elem csv
        with open('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/faults_elem_properties.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # write the header
            # writer.writerow(header)

            # write multiple rows
            writer.writerows(arr_elem_properties)
    elif mesh_type == "leftpostkink":
        #save elem ids
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_mainpath_left.txt',
                    ind_elem_mainpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_branchpath_left.txt',
                    ind_elem_branchpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_daughterpath_left.txt',
                    ind_elem_daughterpath,
                    fmt='%i',
                    newline=" ")

        #save subdomain id
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_mainpath_left.txt',
                    ind_subdomainid_mainpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_branchpath_left.txt',
                    ind_subdomainid_branchpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_daughterpath_left.txt',
                    ind_subdomainid_daughterpath,
                    fmt='%i',
                    newline=" ")

        #save property elem csv
        with open('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/faults_elem_properties_left.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # write the header
            # writer.writerow(header)

            # write multiple rows
            writer.writerows(arr_elem_properties)
    elif mesh_type == "rightpostkink":
        #save elem ids
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_mainpath_right.txt',
                    ind_elem_mainpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_branchpath_right.txt',
                    ind_elem_branchpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_elemid_daughterpath_right.txt',
                    ind_elem_daughterpath,
                    fmt='%i',
                    newline=" ")

        #save subdomain id
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_mainpath_right.txt',
                    ind_subdomainid_mainpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_branchpath_right.txt',
                    ind_subdomainid_branchpath,
                    fmt='%i',
                    newline=" ")
        np.savetxt('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/list_subdomainid_daughterpath_right.txt',
                    ind_subdomainid_daughterpath,
                    fmt='%i',
                    newline=" ")

        #save property elem csv
        with open('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/faults_elem_properties_right.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # write the header
            # writer.writerow(header)

            # write multiple rows
            writer.writerows(arr_elem_properties)

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

def read_properties(xlsx_file_read_path,mesh_type):

    df = pd.read_excel(xlsx_file_read_path)

    arr_data = np.array(df)

    if mesh_type == "full":
        arr_data = arr_data[1:,[4,5,7,8,10]]

    return arr_data


if __name__ == "__main__":

    #--------------------------points/lines generation--------------------------#
    # mesh_type : type of mesh (full, leftpostkink, rightpostkink)
    mesh_type = "full"
    #---------------------------------------------------------------------------#

    scale = 100000
    shriftdist_x = 50.992
    shriftdist_y = 13.9
    # inputfilepath = "/Users/andyz/Documents/MOOSE_FARMS/geometry/vertices.csv"
    inputfilepath = "/Users/andyz/Documents/MOOSE_FARMS/geometry/vertices_rotated.csv"
    my_data = np.genfromtxt(inputfilepath, delimiter=',',skip_header=1)
    data_xcoord = my_data[:,0]
    data_ycoord = my_data[:,1]
    angle = -30

    #rotate the cooridnates
    coordxr, coordyr = helper_rotate(data_xcoord,data_ycoord,angle)

    #physical coordinates
    coordxp, coordyp = helper_physical(coordxr,coordyr,shriftdist_x,shriftdist_y,scale)
    
    #move coordinates
    coordxm, coordym = helper_movecoords(coordxp,coordyp)

    #
    dictlines = {1 : [0,0,-33921.924,-21196.771],
                 2 : [0,0,-30000,0],
                 3 : [-30000,0,-102600,-52061.46],
                 4 : [-102600,-52061.46,-109440,-52760],
                 5 : [-109440,-52760,-117070,-57430],
                 6 : [-113530.2,-55508,-117920,-61040],
                 7 : [-117920,-61040,-124846,-64924],
                 8 : [-139630,-74330,-143617,-81610],
                 9 : [-124846,-64924,-139630,-74330],
                10 : [-132999,-70101,-135970,-75040],
                11 : [0,0,50000,0],
                12 : [50000,0,79240,-11830],
                13 : [50000,0,63500,-6],
                14 : [63500,-6,146440,-17500],
                15 : [-35885.573,-4219.874,-42517.376,-15541.61],
                16 : [-36743.403,-4834.926,-40690,-10000],
                17 : [-65114.711,-25176.7612,-74230,-26200],
                18 : [-139630,-74330,-155520,-83320],
                19 : [-144570,-77130,-155160,-77170],
                20 : [-124846,-64924,-135780,-66170]}
    
    #change 1: segment 6 node coordinates

    # print(len(dictlines))
    #add piecewise constant line
    helper_addlines(coordxm,coordym,dictlines)

    # get points and connectivity for meshing purpose
    dict_unique_ptrs = helper_getallptrs(dictlines, mesh_type=mesh_type)

    # #add intersection points
    # sizeofdictuniqueptrs = len(dict_unique_ptrs)
    # dict_unique_ptrs[sizeofdictuniqueptrs+1] = [-145267.434, -61222.656]

    #construct connectivity
    if mesh_type == "full":
        dict_lines_connectivity = { 1 : [1,3],
                                    2 : [2,3],
                                    3 : [3,5],
                                    4 : [4,5],
                                    5 : [5,8],
                                    6 : [6,8],
                                    7 : [8,9],
                                    8 : [7,9],
                                    9 : [9,10],
                                    10 : [10,12],
                                    11 : [11,13],
                                    12 : [13,14],
                                    13 : [14,16],
                                    14 : [15,16],
                                    15 : [16,19],
                                    16 : [18,19],
                                    17 : [19,20],
                                    18 : [17,20],
                                    19 : [20,22],
                                    20 : [22,23],
                                    21 : [21,23],
                                    22 : [23,24],
                                    23 : [24,26],
                                    24 : [24,25],
                                    25 : [25,27]}
    
    #change 2: delete segment 12, redo numbering after

    #modify points/lines to lower mu_s
    # dict_unique_ptrs[23] = [-74910,-14350]

    #print/check
    helper_addmeshlines(dict_lines_connectivity,dict_unique_ptrs,coordxm,coordym)

    #reverse the coordinates#
    #helper_reverse(coordxp,coordyp,dict_unique_ptrs,dict_lines_connectivity,scale,shriftdist_x,shriftdist_y,angle)
    exit()
    #--------------------------property calculation--------------------------#
    if mesh_type == "full":
        #assume
        #
        ratio_psts1_psts3 = 1/4                #ratio sigma1/sigma3
        psts3 = -15                            #magnitude of sigma3 MPa
        #
        psts1 = psts3 / ratio_psts1_psts3      #magnitude of sigma1
        #
        nuc_len = 1.6                          #nuculeation length 1.6km
        shearmodulus = 32.04                   #shear modulus GPa
        #
        angleofmaxcomp = 44                    #angle of maximum compression w.r.t planar fault (degree)
        #
        # static_friction_mus = 0.7              #static friction mu_s
        # static_friction_mud = 0.2              #dynamic friction mu_d
        ##Note: here we assume static friction and dynamic friction are the same across simulation
        #this makes the shear-to-normal ratio as the only variable
        # arr_frictionmus = np.ones((len(dict_lines_connectivity),1)) * static_friction_mus
        # arr_frictionmud = np.ones((len(dict_lines_connectivity),1)) * static_friction_mud
        #
        ind_mainfault = np.array([[30,32,33]]) - 1         #index of main fault segment
        ind_splayfault = np.array([[31]]) - 1              #index of splay fault segment
        ind_mainfaultbeyondkink = np.array([i for i in range(1,29+1)]) #index of main fault beyond kink segment
        #
        S_mf = 1.0
        S_mfbyk = 2.0
        S_sf = 0.75
        # arr_Sparam = np.ones((len(dict_lines_connectivity),1)) #assume 1 initially
        # arr_Sparam[ind_mainfault] = S_mf
        # arr_Sparam[ind_splayfault] = S_sf
        # arr_Sparam[ind_mainfaultbeyondkink] = S_mfbyk
        #
        #important segments on main fault beyond kink should be critical for choosing mu_s and mu_d
        # ind_importantsegments_mainfaultbykink = [21,22,23,25,27,29]

        #compute angles
        arr_lines_angles = helper_Hangle(dict_lines_connectivity,dict_unique_ptrs)

        #compute stress
        arr_angleofmaxcompwfault, arr_normalsts, arr_shearsts, arr_sheartonormal = propertycalc_sts(psts1,psts3,arr_lines_angles,angleofmaxcomp)

        #save to csv
        #id angle angle_comp shear_sts normal-sts shearnormalratio
        
        header = ["segment id", "angle", "angle max compression", "initial shear stress", "initial normal stress", "shear to normal ratio"]
        
        arr_data = np.zeros((np.shape(arr_angleofmaxcompwfault)[0],6))

        #id
        arr_data[:,0] = np.linspace(1,np.shape(arr_angleofmaxcompwfault)[0],np.shape(arr_angleofmaxcompwfault)[0])

        #angle
        arr_data[:,1] = arr_lines_angles.reshape(np.shape(arr_lines_angles)[0],)

        #angle_comp
        arr_data[:,2] = arr_angleofmaxcompwfault.reshape(np.shape(arr_angleofmaxcompwfault)[0],)

        #shear stress
        arr_data[:,3] = arr_shearsts.reshape(np.shape(arr_shearsts)[0],)

        #normal stress
        arr_data[:,4] = arr_normalsts.reshape(np.shape(arr_normalsts)[0],)

        #shear to normal ratio
        print(arr_sheartonormal)
        arr_data[:,5] = arr_sheartonormal.reshape(np.shape(arr_sheartonormal)[0],)

        with open('/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/faultsdata.csv', 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)

            # write the header
            writer.writerow(header)

            # write multiple rows
            writer.writerows(arr_data)

        #frictional properties are determined in faultsdata.csv

    #show all figures in the end
    plt.show()

    #--------------------------mesh generation--------------------------#

    #mesh size
    lc = 200

    #number meshing algorithm
    num_algo_meshing = 6

    #file path
    if mesh_type == "full":
        file_path = "/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/complexfaults.msh"

    #generate new mesh
    new_mesh = False
    #add corner points
    sizeofdictuniqueptrs = len(dict_unique_ptrs)
    
    if mesh_type == "full":
        dict_unique_ptrs[sizeofdictuniqueptrs+1] = [-158000,-85000]
        dict_unique_ptrs[sizeofdictuniqueptrs+2] = [148500,-85000]
        dict_unique_ptrs[sizeofdictuniqueptrs+3] = [148500,10000]
        dict_unique_ptrs[sizeofdictuniqueptrs+4] = [-158000,10000]

    #add elems
    sizeofdictlinesconnect = len(dict_lines_connectivity)

    dict_lines_connectivity[sizeofdictlinesconnect+1] = [sizeofdictuniqueptrs+1,sizeofdictuniqueptrs+2]
    dict_lines_connectivity[sizeofdictlinesconnect+2] = [sizeofdictuniqueptrs+2,sizeofdictuniqueptrs+3]
    dict_lines_connectivity[sizeofdictlinesconnect+3] = [sizeofdictuniqueptrs+3,sizeofdictuniqueptrs+4]
    dict_lines_connectivity[sizeofdictlinesconnect+4] = [sizeofdictuniqueptrs+4,sizeofdictuniqueptrs+1]

    #add line loop
    line_loop = [sizeofdictlinesconnect+1,sizeofdictlinesconnect+2,sizeofdictlinesconnect+3,sizeofdictlinesconnect+4]

    if new_mesh == True:

        #gmsh add points
        gmsh_add_points(dict_unique_ptrs,lc,mesh_type=mesh_type)

        #gmsh add elems
        gmsh_add_elems(dict_lines_connectivity,mesh_type=mesh_type)

        #gmsh add surfaces
        gmsh_add_surfaces(line_loop)

        #gmsh embed points/lines
        gmsh_embed_entities(dict_unique_ptrs,dict_lines_connectivity,mesh_type=mesh_type)

        #gmsh add physical groups
        gmsh_add_physical_entities(dict_lines_connectivity,mesh_type=mesh_type)

        #gmsh meshing
        gmsh_meshing(num_algo_meshing, file_path)

    #------------------------Generate mesh subdomain IDs------------------------#

    splayfault_nucleation_shear_sts = -12.348 #MPa

    if mesh_type == "full":
        main_path_segmentids = [1,3,5,7,9,10,11,12,13,15,17,19,20,22,24]
        branch_path_segmentids = [2,4,6,8,14,16,18,21,23,25]
        daughter_path_segmentids = []

    # file_path = "/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/complexfaults.msh"
    
    xlsx_file_read_path = "faultsdata_main_excel.xlsx"
    
    #csv_file_save_path = "/Users/andyz/projects/febe_moose/problems_febe/test_CZM_turkeybranch/meshv8full/complexfaults.msh"
    
    #initial shear stress, initial normal stress, mu_s, mu_d, Dc

    header = ["ini_shear_sts","ini_normal_sts","mu_s","mu_d","Dc"]

    arr_data = read_properties(xlsx_file_read_path,mesh_type=mesh_type)

    print(arr_data)

    meshio_subdomainids(dict_lines_connectivity,
                        dict_unique_ptrs,
                        file_path,
                        main_path_segmentids,
                        branch_path_segmentids,
                        daughter_path_segmentids,
                        arr_data,
                        splayfault_nucleation_shear_sts,
                        header)
