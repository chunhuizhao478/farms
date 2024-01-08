"""
Multifault Generation

Functions

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
from scipy.spatial import Delaunay
import numpy as np

## Generate Boundary Ptrs Coordinates
def gen_coord_1st(list_cornerptrs,
                  case_flag):

    """
    [inputs]
    list_cornerptrs: list of corner ptrs coordinates [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]
    case_flag: classify different boundary points definition
    
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

    well:
      - - -
    |       |
    5   6   7
    |       |
      - - - 
    """
    # initialize array_coord
    numofnodes = int(len(list_cornerptrs) / 2)

    array_coord = np.zeros(( 2 * numofnodes,1)) 

    # store four corner + two main fault ptr coord
    for i in range(int(len(list_cornerptrs)/2)):
        array_coord[2 * i] = list_cornerptrs[2 * i]
        array_coord[2 * i + 1] = list_cornerptrs[2 * i + 1]
    
    array_coord = array_coord.reshape((numofnodes,2))

    return array_coord

## GMSH function to create points given coordinates
def gmsh_add_ptrs_1st(array_coord, 
                      lc1,
                      case_flag):

    """
    [inputs]
    array_coord: coordinates of ptr coordinates ([ptri_x,ptri_y])
    lc1: mesh size specified at point
    case_flag: classify different boundary points definition

    [defines]
    x_i : x coordinate for ptr i
    y_i : y coordinate for ptr i
    p1_x : x coordinate for ptr 1
    p1_y : y coordinate for ptr 1
    p2_x : x coordinate for ptr 2
    p2_y : y coordinate for ptr 2
    list_gmsh_corner_ptr : temporary list for storing gmsh point identities for corner/edge points
    additional_files : tuple of storing any additional structure 

    Note:
    additional_files["circle"] for access
    
    [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]
    
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

   well:
      - - -
    |       |
    5   6   7
    |       |
      - - - 
    """

    #gmsh initialize
    gmsh.initialize()

    #add corner points
    list_gmsh_corner_ptr = []
    for corner_ptr_ind in range(np.shape(array_coord)[0]):
        x_i = array_coord[corner_ptr_ind, 0]
        y_i = array_coord[corner_ptr_ind, 1]
        ptr_i = gmsh.model.occ.addPoint(x_i,y_i,0,lc1)
        # print(ptr_i)
        list_gmsh_corner_ptr.append(ptr_i)
    print("finish adding " + str(np.shape(array_coord)[0]) + " corner points")

    #add to mat_connect_corner
    """
    Follow the note in [defines] section, connectivity is given as follows:
    1-2, 2-3, 3-4, 4-1
    case_flag: 2D-Cluster

    case_flag: 2D-Cluster-Well
    other than boundary, well:
    5-6-7, 7-6-5
    """

    #output additional files
    additional_mats = {}

    if case_flag == "2D-Cluster" or "2D-Cluster-Well":
        mat_connect_corner = np.zeros((4,2),dtype=int)
        mat_connect_corner[0,0] = list_gmsh_corner_ptr[0]; mat_connect_corner[0,1] = list_gmsh_corner_ptr[1]
        mat_connect_corner[1,0] = list_gmsh_corner_ptr[1]; mat_connect_corner[1,1] = list_gmsh_corner_ptr[2]
        mat_connect_corner[2,0] = list_gmsh_corner_ptr[2]; mat_connect_corner[2,1] = list_gmsh_corner_ptr[3]
        mat_connect_corner[3,0] = list_gmsh_corner_ptr[3]; mat_connect_corner[3,1] = list_gmsh_corner_ptr[0]

    #add circle arc 
    if case_flag == "2D-Cluster-Well":
        mat_connect_circlearc = np.zeros((2,3),dtype=int)   
        mat_connect_circlearc[0,0] = list_gmsh_corner_ptr[4]; mat_connect_circlearc[0,1] = list_gmsh_corner_ptr[5]; mat_connect_circlearc[0,2] = list_gmsh_corner_ptr[6]
        mat_connect_circlearc[1,0] = list_gmsh_corner_ptr[6]; mat_connect_circlearc[1,1] = list_gmsh_corner_ptr[5]; mat_connect_circlearc[1,2] = list_gmsh_corner_ptr[4]

        additional_mats["circle"] = mat_connect_circlearc

    return mat_connect_corner, additional_mats

# GMSH function to create elements given connectivity
def gmsh_add_elems_1st(mat_connect_corner,
                       additional_mats,
                       case_flag):

    """
    [inputs]
    mat_connect_corner : connectivity matrix for corner points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
    additional_mats : tuple of storing any additional structure 
    Note:
    additional_mats["circle"] for access

    [defines]
    ptr_i : GMSH ptr i (integar)
    ptr_j : GMSH ptr j (integar)
    arr_elem_corner : array stores GMSH element id for edge elements   [elem_k]
    mat_surface_lineloop : matrix for storing surface connectivity [elem1; elem2; elem3; elem4] : elem1 -> elem2 -> elem3 -> elem4 -> elem1
    arr_elem_additional: array stores GMSH element id (other than line)  [elem_k], arr_elem_additional["circle"] = a list of element id
    mat_surface_additional_lineloop: additional surface loop (other than line)

    case_flag: 2D-Cluster
    1-2, 2-3, 3-4, 4-1
    4 - - -[3] - - 3
    |              |
    |              |
    |              |
   [4]            [2]
    |              |
    |              |
    |              |
    1 - - [1]- - - 2

    case_flag: 2D-Cluster-Well
    1-2, 2-3, 3-4, 4-1
    4 - - -[3] - - 3
    |              |
    |              |
    |              |
   [4]      o     [2]
    |              |
    |              |
    |              |
    1 - - [1]- - - 2

    well:
      - - -
    |       |
   [5] [6] [7]
    |       |
      - - - 

    """

    #initialize arr_elem_corner
    arr_elem_corner = np.zeros((np.shape(mat_connect_corner)[0],1),dtype=int)

    #additional elems type
    arr_elem_additional = {}

    #add edge elements 
    for corner_elem_ind in range(np.shape(mat_connect_corner)[0]):
        ptr_i = mat_connect_corner[corner_elem_ind,0]
        ptr_j = mat_connect_corner[corner_elem_ind,1]
        elem_k = gmsh.model.occ.addLine(ptr_i,ptr_j)
        arr_elem_corner[corner_elem_ind] = elem_k
    print("finish adding " + str(np.shape(mat_connect_corner)[0]) + " edge elements")

    #add circle arcs
    if case_flag == "2D-Cluster-Well":
        arr_elem_additional["circle"] = np.zeros((np.shape(additional_mats["circle"])[0],1),dtype=int)
        for arc_elem_ind in range(np.shape(additional_mats["circle"])[0]):
            ptr_start = additional_mats["circle"][arc_elem_ind,0]
            ptr_center = additional_mats["circle"][arc_elem_ind,1]
            ptr_end = additional_mats["circle"][arc_elem_ind,2]
            circlearc_k = gmsh.model.occ.addCircleArc(ptr_start,ptr_center,ptr_end)
            arr_elem_additional["circle"][arc_elem_ind] = circlearc_k
    
    #additional surface loop (other than line)
    mat_surface_additional_lineloop = {}

    #add mat_surface_lineloop
    if case_flag == "2D-Cluster" or "2D-Cluster-Well" :
        mat_surface_lineloop = np.zeros((4,1),dtype=int)
        ## only one surface
        mat_surface_lineloop[0,0] = arr_elem_corner[0];
        mat_surface_lineloop[1,0] = arr_elem_corner[1];
        mat_surface_lineloop[2,0] = arr_elem_corner[2];
        mat_surface_lineloop[3,0] = arr_elem_corner[3];
    
    if case_flag == "2D-Cluster-Well":
        mat_surface_additional_lineloop["circle"] = np.zeros((2,1),dtype=int)
        ## only one circle
        mat_surface_additional_lineloop["circle"][0,0] = arr_elem_additional["circle"][0]
        mat_surface_additional_lineloop["circle"][1,0] = arr_elem_additional["circle"][1]

    return mat_surface_lineloop, arr_elem_corner, mat_surface_additional_lineloop, arr_elem_additional

# GMSH function to create surface given line loop
def gmsh_add_surfaces_1st(mat_surface_lineloop,
                          mat_surface_additional_lineloop,
                          case_flag):

    """
    [inputs]
    mat_surface_lineloop : matrix for storing surface connectivity [elem1; elem2; elem3; elem4] : elem1 -> elem2 -> elem3 -> elem4 -> elem1
    mat_surface_additional_lineloop : tuple for storing surface connectivity : mat_surface_additional_lineloop["circle"]
    case_flag : "2D-Cluster",  "2D-Cluster-Well"
    [defines]
    loopi : line loop i
    surfi : surface i

    [returns]
    arr_surface :  array stores GMSH id of surface
    arr_additional_surface : array stores GMSH id of non-straight-line-loop surface
    """

    #get number of surfaces
    num_of_surfaces = np.shape(mat_surface_lineloop)[1]

    #initialize arr_surface
    arr_surface = np.zeros((num_of_surfaces,1),dtype=int)

    #additional surface
    arr_additional_surface = {}

    #add surfaces
    for surface_ind in range(num_of_surfaces):
        loopi = gmsh.model.occ.addCurveLoop([mat_surface_lineloop[0,surface_ind],
                                             mat_surface_lineloop[1,surface_ind],       
                                             mat_surface_lineloop[2,surface_ind],
                                             mat_surface_lineloop[3,surface_ind]])
        if case_flag == "2D-Cluster":
            surfi = gmsh.model.occ.addPlaneSurface([loopi])
            arr_surface[surface_ind] = surfi
    
    #add circle surfaces
    if case_flag == "2D-Cluster-Well":
        arr_additional_surface["circle"] = np.zeros((np.shape(mat_surface_additional_lineloop["circle"])[1],1),dtype=int)
        for circle_surface_ind in range(np.shape(mat_surface_additional_lineloop["circle"])[1]):
            #only one circle
            circle_loopi = gmsh.model.occ.addCurveLoop([mat_surface_additional_lineloop["circle"][0,circle_surface_ind],
                                                        mat_surface_additional_lineloop["circle"][1,circle_surface_ind]],100)
            circle_surface = gmsh.model.occ.addPlaneSurface([loopi,circle_loopi]) #one surface, one circle, hardcode
            arr_additional_surface["circle"][0] = circle_surface

    #synchronize
    # gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    return arr_surface, arr_additional_surface

# GMSH function to create physical entities
def gmsh_add_physical_entities_1st(arr_elem_corner, 
                                   arr_elem_additional,
                                   arr_surface,
                                   arr_additional_surface,
                                   case_flag):

    """
    [inputs]
    arr_elem_corner : array stores GMSH element id for edge elements   [elem_k]
    arr_surface : array stores GMSH id of surface
    arr_elem_additional : array stores GMSH element id for edge elements   [elem_k]
    arr_additional_surface : array stores GMSH id of non-line-loop surface

    case_flag: 2D-Cluster, 2D-Cluster-Well
    1-2, 2-3, 3-4, 4-1
    4 - - -[3] - - 3
    |              |
    |              |
    |              |
   [4]            [2]
    |              |
    |              |
    |              |
    1 - - [1]- - - 2

    physical curves:
    top: 3
    bottom: 1
    left: 4
    right: 2

    additional:
    case_flag: 2D-Cluster-Well

     well:
      - - -
    |       |
   [5] [6] [7]
    |       |
      - - - 

    [returns]
    mat_egelem_phyname: record edge element's physical name [elem_gmsh_id, elem_phy_name]
    mat_surf_phyname: record surface's physical name [surf_gmsh_id, surf_phy_name] 
    """

    #initialize mat_elem_phyname, mat_surf_phyname
    mat_egelem_phyname = {}
    mat_surf_phyname = {}

    mat_egelem_phyname_additional = {} #note the key/value pair is different from above
    mat_surface_phyname_additional = {} #note the key/value pair is different from above

    #add boundary physical group
    if case_flag == "2D-Cluster" or "2D-Cluster-Well":
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[2,0]],name="top")
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[0,0]],name="bottom")
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[3,0]],name="left")
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[1,0]],name="right")

        #save in mat_egelem_phyname
        mat_egelem_phyname[arr_elem_corner[2,0]] = "top"
        mat_egelem_phyname[arr_elem_corner[0,0]] = "bottom"
        mat_egelem_phyname[arr_elem_corner[3,0]] = "left"
        mat_egelem_phyname[arr_elem_corner[1,0]] = "right"
    
    #add circlearc physical group
    if case_flag == "2D-Cluster-Well":
        gmsh.model.addPhysicalGroup(1,[arr_elem_additional["circle"][0,0]],name="well")
        gmsh.model.addPhysicalGroup(1,[arr_elem_additional["circle"][1,0]],name="well")
        mat_egelem_phyname_additional["well"] = [arr_elem_additional["circle"][0,0],arr_elem_additional["circle"][1,0]]

    #add surface physical group
    if case_flag == "2D-Cluster":
        for surface_id in range(np.shape(arr_surface)[0]):
            gmsh.model.addPhysicalGroup(2,[arr_surface[surface_id,0]],name="surf" + str(arr_surface[surface_id,0]))
            #sabe in mat_surf_phyname
            mat_surf_phyname["surf"] = "surf" + str(arr_surface[surface_id,0])
    elif case_flag == "2D-Cluster-Well":
        for additional_surface_id in range(np.shape(arr_additional_surface["circle"])[0]):
            gmsh.model.addPhysicalGroup(2,[arr_additional_surface["circle"][additional_surface_id,0]],name="surfwell" + str(arr_additional_surface["circle"][additional_surface_id,0]))
            #sabe in mat_surf_phyname
            mat_surface_phyname_additional["surfwell"] = "surfwell" + str(arr_additional_surface["circle"][additional_surface_id,0])

    #synchronize
    gmsh.model.occ.synchronize()

    return mat_egelem_phyname, mat_surf_phyname, mat_egelem_phyname_additional, mat_surface_phyname_additional

# GMSH function to do meshing using provided algorithm 
def gmsh_meshing(num_alg_meshing,file_path1):

    """
    [inputs]
    num_alg_meshing: mesh algorithm number
    file_path: path for msh file to be saved

    list of mesh algorithms:

    Mesh.Algorithm
    1 : MeshAdpt
    2 : Automatic
    3 : Initial mesh only
    5 : Delaunay
    6 : Frontal-Delaunay
    7 : BAMG
    8 : Frontal-Delaunay for Quads
    9 : Packing of Parallelograms
    11 : Quasi-Structured Quad

    Mesh.RecombinationAlgorithm
    0 : simple
    1 : blossom
    2 : simple full-quad
    3 : blossom full-quad
    """
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.Algorithm", num_alg_meshing)
    
    #generate mesh
    gmsh.model.mesh.generate(dim=2)
    gmsh.fltk.run()
    gmsh.write(file_path1)
    gmsh.finalize()

# Save Network Ptrs/Connectivity
def meshio_save_network_info(file_path1,
                             array_coord,
                             dist_verticalbc,
                             dist_lateralbc,
                             case_flag,
                             surfnames,
                             surfnames_additional,
                             radius):

    """
    [inputs]
    file_path1: path for msh file to be saved
    array_coord: coordinates of ptr coordinates ([ptri_x,ptri_y])
    dist_verticalbc: distance in y direction of the network each side
    dist_lateralbc: distance in x diraction of the network each side
    
    if case_flag == 2D-Cluster
    array_coord: list of corner ptrs coordinates [x1, y1, x2, y2, x3, y3, x4, y4]
    4 - - - - - 3
    |           |
    |           |
    |           |
    1 - - - - - 2

    [defines]
    mat_surf1_connect: Define Larger mat for store connectivity [ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y]
    mat_surf1_triacoordall: Define Larger mat for store triangle three coordinates [ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y, ptr_k_x, ptr_k_y]
    mat_surf_ptr_all: Define Larger mat for store branch coordinate
    (Note: end ptrs are saved!)
    """

    m = meshio.read(file_path1)
    
    #Get Coord for all points
    coord_all = m.points

    #Get TRIA3 Connectivity
    tria_elem_connect = m.cells_dict['triangle']

    #Bounds
    if case_flag == "2D-Cluster":
        bound_min_x = array_coord[0,0] #x1
        bound_max_x = array_coord[1,1] #x2
        bound_min_y = array_coord[0,1] #y1
        bound_max_y = array_coord[2,1] #y3

    ####----------------------------SURFACE----------------------------####
    #cases only single surface
    if case_flag == "2D-Cluster" or "2D-Cluster-Well":

        #Get Tri3 Index
        if case_flag == "2D-Cluster":
            surf_tria_ind = m.cell_sets_dict[surfnames['surf']]['triangle'] #hardcode
        elif case_flag == "2D-Cluster-Well":
            surf_tria_ind = m.cell_sets_dict[surfnames_additional['surfwell']]['triangle']

        #Define Larger mat for store connectivity [ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y]
        mat_surf1_connect = np.zeros((np.shape(surf_tria_ind)[0] * 3,4))
        mat_surf1_connect_count = 0

        #Define Larger mat for store triangle three coordinates
        #[ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y, ptr_k_x, ptr_k_y]
        mat_surf1_triacoordall = np.zeros((np.shape(surf_tria_ind)[0] * 3,6))
        mat_surf1_tria_count = 0

        #Define line coord list
        arr_line_coord_x = np.zeros((3,2))
        arr_line_coord_y = np.zeros((3,2))

        #Loop over tria
        for tria_ind in range(np.shape(surf_tria_ind)[0]):

            #get current tria ind
            tria_i = surf_tria_ind[tria_ind]

            #get connectivity for current element
            tria_connect_i = tria_elem_connect[tria_i,:]

            #get coordinates
            coord_x_i = coord_all[tria_connect_i, 0]
            coord_y_i = coord_all[tria_connect_i, 1]

            #construct line connect
            arr_line_coord_x[0,0] = coord_x_i[0]; arr_line_coord_x[0,1] = coord_x_i[1]
            arr_line_coord_x[1,0] = coord_x_i[1]; arr_line_coord_x[1,1] = coord_x_i[2]
            arr_line_coord_x[2,0] = coord_x_i[2]; arr_line_coord_x[2,1] = coord_x_i[0]

            arr_line_coord_y[0,0] = coord_y_i[0]; arr_line_coord_y[0,1] = coord_y_i[1]
            arr_line_coord_y[1,0] = coord_y_i[1]; arr_line_coord_y[1,1] = coord_y_i[2]
            arr_line_coord_y[2,0] = coord_y_i[2]; arr_line_coord_y[2,1] = coord_y_i[0]

            #check1: within the network bounds
            if (np.max(arr_line_coord_y) < dist_verticalbc and np.min(arr_line_coord_y) > -dist_verticalbc and np.min(arr_line_coord_x) > -1 * dist_lateralbc and np.max(arr_line_coord_x) < 1 * dist_lateralbc ):

                #loop over all three line elements
                for line_ind in range(np.shape(arr_line_coord_x)[0]):

                    #get current line coord
                    arr_line_coord_x_current = arr_line_coord_x[line_ind,:]
                    arr_line_coord_y_current = arr_line_coord_y[line_ind,:]
                        
                    #save in mat_surf1_connect
                    mat_surf1_connect[mat_surf1_connect_count, 0] = arr_line_coord_x_current[0]
                    mat_surf1_connect[mat_surf1_connect_count, 1] = arr_line_coord_y_current[0]
                    mat_surf1_connect[mat_surf1_connect_count, 2] = arr_line_coord_x_current[1]
                    mat_surf1_connect[mat_surf1_connect_count, 3] = arr_line_coord_y_current[1]

                    #update count
                    mat_surf1_connect_count += 1
                
                #hardcode

                #save this element all coordinates
                mat_surf1_triacoordall[mat_surf1_tria_count, 0] = coord_x_i[0]
                mat_surf1_triacoordall[mat_surf1_tria_count, 1] = coord_y_i[0]
                mat_surf1_triacoordall[mat_surf1_tria_count, 2] = coord_x_i[1]
                mat_surf1_triacoordall[mat_surf1_tria_count, 3] = coord_y_i[1]
                mat_surf1_triacoordall[mat_surf1_tria_count, 4] = coord_x_i[2]
                mat_surf1_triacoordall[mat_surf1_tria_count, 5] = coord_y_i[2]

                #update count
                mat_surf1_tria_count += 1

            else:
                continue

        #eliminate zero entries for mat
        mat_surf1_connect = mat_surf1_connect[:mat_surf1_connect_count,:]

        mat_surf1_triacoordall = mat_surf1_triacoordall[:mat_surf1_tria_count,:]

    ####------------------------------------------------------------------####

    if case_flag == "2D-Cluster" or "2D-Cluster-Well":
        #Attach two mat_surf_connect
        mat_surf_connect_all = mat_surf1_connect

        #Attach two mat_surf_triacoordall
        mat_surf_triacoordall = mat_surf1_triacoordall

    #Extend the mat_surf_triacoordall into mat_surf_trialinesegmentall 
    #[ptr_i_x ptr_i_y ptr_j_x ptr_j_y ptr_j_x ptr_j_y ptr_k_x ptr_k_y ptr_k_x ptr_k_y ptr_i_x ptr_i_y]
    #such that line segement i-j j-k k-i
    mat_surf_trialinesegmentall = np.zeros((np.shape(mat_surf_triacoordall)[0],12))
    mat_surf_trialinesegmentall[:,0] = mat_surf_triacoordall[:,0]
    mat_surf_trialinesegmentall[:,1] = mat_surf_triacoordall[:,1]
    mat_surf_trialinesegmentall[:,2] = mat_surf_triacoordall[:,2]
    mat_surf_trialinesegmentall[:,3] = mat_surf_triacoordall[:,3]
    mat_surf_trialinesegmentall[:,4] = mat_surf_triacoordall[:,2]
    mat_surf_trialinesegmentall[:,5] = mat_surf_triacoordall[:,3]
    mat_surf_trialinesegmentall[:,6] = mat_surf_triacoordall[:,4]
    mat_surf_trialinesegmentall[:,7] = mat_surf_triacoordall[:,5]
    mat_surf_trialinesegmentall[:,8] = mat_surf_triacoordall[:,4]
    mat_surf_trialinesegmentall[:,9] = mat_surf_triacoordall[:,5]
    mat_surf_trialinesegmentall[:,10] = mat_surf_triacoordall[:,0]
    mat_surf_trialinesegmentall[:,11] = mat_surf_triacoordall[:,1]

    #Reshape mat_suf_connect_all to (shape_x,2)
    mat_surf_connect_all_reshape = np.reshape(mat_surf_connect_all,(np.shape(mat_surf_connect_all)[0] * 2,2))

    #Unique mat_surf_connect_all_reshape
    #Define Larger mat for store branch coordinate
    mat_surf_ptr_all = np.unique(mat_surf_connect_all_reshape,axis=0)

    # Initialize mat_connect_additional [ptri_x ptri_y ptrj_x ptrj_y] line segment i -> j
    mat_connect_additional = {}
    mat_ptrs_additional = {}

    ## Find Approximate line segments surrounding the circle ##
    if case_flag == "2D-Cluster-Well":
        list_circlebc = [] 
        # Loop Over point coordinates
        for ptr_ind in range(np.shape(coord_all)[0]):
            ptr_i_x = coord_all[ptr_ind,0]
            ptr_i_y = coord_all[ptr_ind,1]
            if ptr_i_x * ptr_i_x + ptr_i_y * ptr_i_y - radius * radius < 1e-1:
                list_circlebc.append(ptr_i_x)
                list_circlebc.append(ptr_i_y)
        #asarray
        arr_circlebe = np.array(list_circlebc)
        arr_circlebe = np.reshape(arr_circlebe,[int(np.shape(arr_circlebe)[0]/2),2])
        arr_circlebe_x = arr_circlebe[:,0]
        arr_circlebe_y = arr_circlebe[:,1]
        arr_origin_x = np.zeros(np.shape(arr_circlebe_x))
        arr_origin_y = np.zeros(np.shape(arr_circlebe_y))
        #determine sequence #hardcode origin (0,0)
        angles = helper_angle_between_points(arr_circlebe_x, arr_circlebe_y, arr_origin_x, arr_origin_y)
        #sort 
        sorted_ind = np.argsort(angles)
        arr_circlebe_x_sorted = arr_circlebe_x[sorted_ind]
        arr_circlebe_y_sorted = arr_circlebe_y[sorted_ind]
        #add initial one at end
        arr_circlebe_x_sorted = np.append(arr_circlebe_x_sorted,arr_circlebe_x_sorted[0])
        arr_circlebe_y_sorted = np.append(arr_circlebe_y_sorted,arr_circlebe_y_sorted[0])
        #save results
        mat_connect_additional["pseudo_circle"] = np.zeros((np.shape(arr_circlebe)[0],4))
        for line_i in range(np.shape(arr_circlebe)[0]):
            mat_connect_additional["pseudo_circle"][line_i,0] = arr_circlebe_x_sorted[line_i]
            mat_connect_additional["pseudo_circle"][line_i,1] = arr_circlebe_y_sorted[line_i]
            mat_connect_additional["pseudo_circle"][line_i,2] = arr_circlebe_x_sorted[line_i+1]
            mat_connect_additional["pseudo_circle"][line_i,3] = arr_circlebe_y_sorted[line_i+1]

        #save points
        mat_ptrs_additional["pseudo_circle"] = arr_circlebe

        #check and remove elements that will be arc (not line)
        # mat_surf_triacoordall[:,0:2]

        # ## for outer faults boundary ##
        # # Computing the alpha shape
        # #hardcode
        # points = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv6_cluster/meshv6/faultpoints.txt")
        # points = np.array(points).reshape([int(len(points)/2),2])

        # edges = helper_alpha_shape(points, alpha=radius, only_outer=True)

        # print(edges)
        # print(len(edges))

        # # Plotting the output
        # plt.figure()
        # plt.axis('equal')
        # plt.plot(points[:, 0], points[:, 1], '.')
        # for i, j in edges:
        #     plt.plot(points[[i, j], 0], points[[i, j], 1])
        # plt.show()

        # #get unique points
        # ptrs_outerbc_list = []
        # for i, j in edges:
        #     ptrs_outerbc_list.append(i)
        #     ptrs_outerbc_list.append(j)
        # arr_ptrs_outerbc = np.unique(ptrs_outerbc_list)

        # print(arr_ptrs_outerbc)

        # # get points coordinates
        # points_outerbc = mat_surf_ptr_all[arr_ptrs_outerbc,:] #hardcode this is still mat_surf_ptr_all

        # # zero as initial points
        # arr_origin_x = np.zeros(np.shape(points_outerbc)[0])
        # arr_origin_y = np.zeros(np.shape(points_outerbc)[0])
    
        # #determine sequence #hardcode origin (0,0)
        # angles = helper_angle_between_points(points_outerbc[:,0], points_outerbc[:,1], arr_origin_x, arr_origin_y)
    
        # #sort 
        # sorted_ind = np.argsort(angles)
        # outerbc_x_sorted = points_outerbc[:,0][sorted_ind]
        # outerbc_y_sorted = points_outerbc[:,1][sorted_ind]
    
        # #add initial one at end
        # outerbc_x_sorted = np.append(outerbc_x_sorted,outerbc_x_sorted[0])
        # outerbc_y_sorted = np.append(outerbc_y_sorted,outerbc_y_sorted[0])

        # # Generate line segement matrix [ptri_x ptri_y ptrj_x ptrj_y] line segment i -> j
        # mat_connect_additional["fault_outerbc"] = np.zeros((len(edges),4))
        # for line_i in range(np.shape(points_outerbc)[0]):
        #     mat_connect_additional["fault_outerbc"][line_i,0] = outerbc_x_sorted[line_i]
        #     mat_connect_additional["fault_outerbc"][line_i,1] = outerbc_y_sorted[line_i]
        #     mat_connect_additional["fault_outerbc"][line_i,2] = outerbc_x_sorted[line_i+1]
        #     mat_connect_additional["fault_outerbc"][line_i,3] = outerbc_y_sorted[line_i+1]
        
        # mat_ptrs_additional["fault_outerbc"] = points_outerbc
        
    #plot
    plt.figure(figsize = (40,4))
    for i in range(np.shape(mat_surf_connect_all)[0]):
        plt.plot([mat_surf_connect_all[i,0],mat_surf_connect_all[i,2]],[mat_surf_connect_all[i,1],mat_surf_connect_all[i,3]],'k-')
    for i in range(np.shape(mat_surf_ptr_all)[0]):
        plt.plot(mat_surf_ptr_all[i,0],mat_surf_ptr_all[i,1],'b*')
    for k in range(np.shape(mat_connect_additional["pseudo_circle"])[0]):
        plt.plot([mat_connect_additional["pseudo_circle"][k,0],mat_connect_additional["pseudo_circle"][k,2]],[mat_connect_additional["pseudo_circle"][k,1],mat_connect_additional["pseudo_circle"][k,3]],'r-')
    plt.show()

    # exit()

    return mat_surf_ptr_all, mat_surf_connect_all, mat_surf_triacoordall, mat_surf_trialinesegmentall, mat_ptrs_additional, mat_connect_additional

def helper_angle_between_points(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    angles = np.degrees(np.arctan2(dy, dx)) % 360
    return angles

#help feed in 2nd circle with smaller radius
def helper_modify_array_coord(array_coord,
                              list_circleptrs_2nd,
                              case_flag):
    
    if case_flag == "2D-Cluster-Well":
        array_coord[4,0] = list_circleptrs_2nd[0];  array_coord[4,1] = list_circleptrs_2nd[1]; 
        array_coord[5,0] = list_circleptrs_2nd[2];  array_coord[5,1] = list_circleptrs_2nd[3];
        array_coord[6,0] = list_circleptrs_2nd[4];  array_coord[6,1] = list_circleptrs_2nd[5];

    return array_coord

# Generate random fault end ptr coordinates
def gmsh_add_points_2nd(array_coord, 
                        mat_surf_ptr,
                        mat_surf_connect,
                        mat_surf_trialinesegmentall,
                        mat_connect_additional,
                        lc2,
                        case_flag):

    """
    [inputs]
    mat_surface_lineloop : matrix for storing surface connectivity [elem1 elem2 elem3 elem4] : elem1 -> elem2 -> elem3 -> elem4 -> elem1
    mat_surf_ptr_all : matrix for storing network point coordinates [ptr_x, ptr_y]
    mat_surf_ptr_all_connect : matrix for storing "coordinate values" for line elements [ptr_i_x, ptr_i_y, ptr_j_x, ptr_y_y] ptr_i connects to ptr_j
    mat_surf_trialinesegmentall : matrix for storing line segment coordinate
    mat_connect_additional : matrix for storing non-square geometry connectivity
    - "pseudocircle" [ptri_x ptri_y ptrj_x ptrj_y] line segment i -> j
    case_flag : case specific

    [ptr_i_x ptr_i_y ptr_j_x ptr_j_y ptr_j_x ptr_j_y ptr_k_x ptr_k_y ptr_k_x ptr_k_y ptr_i_x ptr_i_y]
    such that line segement i-j j-k k-i

    lc2: mesh specified at point
   
    [defines]
    list_gmsh_corner_ptr : temporary list for storing gmsh point identities for corner/edge points
    
    [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6]

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

   well:
      - - -
    |       |
    5   6   7
    |       |
      - - - 

    [returns]
    mat_connect_corner : connectivity matrix for corner points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
    mat_connect_mf : connectivity matrix for main fault points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
    mat_connect_surf : connectivity matrix for branch points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
    additional_mats : tuple of storing any additional structure 
    Note:
    additional_mats["circle"] for access
    """

    #gmsh initialize
    gmsh.initialize()
    gmsh.clear()

    #add corner points
    if case_flag == "2D-Cluster" or "2D-Cluster-Well":
        list_gmsh_corner_ptr = []
        for corner_ptr_ind in range(np.shape(array_coord)[0]): 
            x_i = array_coord[corner_ptr_ind, 0]
            y_i = array_coord[corner_ptr_ind, 1]
            ptr_i = gmsh.model.occ.addPoint(x_i,y_i,0,lc2)
            list_gmsh_corner_ptr.append(ptr_i)
        print("finish adding " + str(np.shape(array_coord)[0]) + " corner points")

    #add to mat_connect_corner
    """
    Follow the note in [defines] section, connectivity is given as follows:
    1-2, 2-3, 3-4, 4-1
    case_flag: 2D-Cluster

    case_flag: 2D-Cluster-Well
    other than boundary, well:
    5-6-7, 7-6-5
    """

    #output additional files
    additional_mats = {}

    if case_flag == "2D-Cluster" or "2D-Cluster-Well":
        mat_connect_corner = np.zeros((4,2),dtype=int)
        mat_connect_corner[0,0] = list_gmsh_corner_ptr[0]; mat_connect_corner[0,1] = list_gmsh_corner_ptr[1]
        mat_connect_corner[1,0] = list_gmsh_corner_ptr[1]; mat_connect_corner[1,1] = list_gmsh_corner_ptr[2]
        mat_connect_corner[2,0] = list_gmsh_corner_ptr[2]; mat_connect_corner[2,1] = list_gmsh_corner_ptr[3]
        mat_connect_corner[3,0] = list_gmsh_corner_ptr[3]; mat_connect_corner[3,1] = list_gmsh_corner_ptr[0]

    #add circle arc 
    if case_flag == "2D-Cluster-Well":
        mat_connect_circlearc = np.zeros((2,3),dtype=int)   
        mat_connect_circlearc[0,0] = list_gmsh_corner_ptr[4]; mat_connect_circlearc[0,1] = list_gmsh_corner_ptr[5]; mat_connect_circlearc[0,2] = list_gmsh_corner_ptr[6]
        mat_connect_circlearc[1,0] = list_gmsh_corner_ptr[6]; mat_connect_circlearc[1,1] = list_gmsh_corner_ptr[5]; mat_connect_circlearc[1,2] = list_gmsh_corner_ptr[4]

        additional_mats["circle"] = mat_connect_circlearc

        #initialize additional_mats["pseudo_circle_gmshptrs_connect"] [ptr_i, ptr_j] define line segment i -> j
        additional_mats["pseudo_circle_gmshptrs_connect"] = np.zeros((np.shape(mat_connect_additional["pseudo_circle"])[0],2))
        # additional_mats["faults_outerbc_gmshptrs_connect"] = np.zeros((np.shape(mat_connect_additional["fault_outerbc"])[0],2))

    #initialize connect matrix
    mat_connect_surf = np.zeros((np.shape(mat_surf_connect)[0],2),dtype=int)

    #initialize connect matrix [ptr_i ptr_j ptr_j ptr_k ptr_k ptr_i] such that i-j j-k k-i
    mat_connect_lineendptrs = np.zeros((np.shape(mat_surf_trialinesegmentall)[0],6),dtype=int)

    #add network fault points
    # np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv6_cluster/meshv6/faultpoints",mat_surf_ptr,fmt='%i',newline=" ")
    # exit()

    #loop over all main fault points
    for surf_point_ind in range(np.shape(mat_surf_ptr)[0]):
        
        # get coordinate for current point
        surf_x_i = mat_surf_ptr[surf_point_ind, 0]
        surf_y_i = mat_surf_ptr[surf_point_ind, 1]
        
        # add gmsh entity
        ptr_i = gmsh.model.occ.addPoint(surf_x_i,surf_y_i,0,lc2)

        # replace all the location in mat_surf_connect with same coordinate by ptr_i in mat_connect_surf
        # first ptr
        suitable_pos_ind_0 = np.where(surf_x_i == mat_surf_connect[:,0])[0]
        suitable_pos_ind_1 = np.where(surf_y_i == mat_surf_connect[:,1])[0]

        # find common flag
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_0,suitable_pos_ind_1)

        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_surf[flag_i, 0] = ptr_i

        # second ptr
        suitable_pos_ind_2 = np.where(surf_x_i == mat_surf_connect[:,2])[0]
        suitable_pos_ind_3 = np.where(surf_y_i == mat_surf_connect[:,3])[0]

        # find common flag
        common_flag_ptr_j = np.intersect1d(suitable_pos_ind_2,suitable_pos_ind_3)

        if common_flag_ptr_j.size > 0:
            for flag_ind_j in range(len(common_flag_ptr_j)):
                flag_j = common_flag_ptr_j[flag_ind_j]
                mat_connect_surf[flag_j,1] = ptr_i
        
        # mat_connect_lineendptrs
        suitable_pos_ind_1 = np.where(surf_x_i == mat_surf_trialinesegmentall[:,0])[0]
        suitable_pos_ind_2 = np.where(surf_y_i == mat_surf_trialinesegmentall[:,1])[0]
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_1,suitable_pos_ind_2)
        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_lineendptrs[flag_i, 0] = ptr_i
        suitable_pos_ind_3 = np.where(surf_x_i == mat_surf_trialinesegmentall[:,2])[0]
        suitable_pos_ind_4 = np.where(surf_y_i == mat_surf_trialinesegmentall[:,3])[0]
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_3,suitable_pos_ind_4)
        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_lineendptrs[flag_i, 1] = ptr_i
        suitable_pos_ind_5 = np.where(surf_x_i == mat_surf_trialinesegmentall[:,4])[0]
        suitable_pos_ind_6 = np.where(surf_y_i == mat_surf_trialinesegmentall[:,5])[0]
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_5,suitable_pos_ind_6)
        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_lineendptrs[flag_i, 2] = ptr_i
        suitable_pos_ind_7 = np.where(surf_x_i == mat_surf_trialinesegmentall[:,6])[0]
        suitable_pos_ind_8 = np.where(surf_y_i == mat_surf_trialinesegmentall[:,7])[0]
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_7,suitable_pos_ind_8)
        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_lineendptrs[flag_i, 3] = ptr_i
        suitable_pos_ind_9 = np.where(surf_x_i == mat_surf_trialinesegmentall[:,8])[0]
        suitable_pos_ind_10 = np.where(surf_y_i == mat_surf_trialinesegmentall[:,9])[0]
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_9,suitable_pos_ind_10)
        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_lineendptrs[flag_i, 4] = ptr_i
        suitable_pos_ind_11 = np.where(surf_x_i == mat_surf_trialinesegmentall[:,10])[0]
        suitable_pos_ind_12 = np.where(surf_y_i == mat_surf_trialinesegmentall[:,11])[0]
        common_flag_ptr_i = np.intersect1d(suitable_pos_ind_11,suitable_pos_ind_12)
        if common_flag_ptr_i.size > 0:
            for flag_ind_i in range(len(common_flag_ptr_i)):
                flag_i = common_flag_ptr_i[flag_ind_i]
                mat_connect_lineendptrs[flag_i, 5] = ptr_i
        
        #additional segments check
        if case_flag == "2D-Cluster-Well":
            #find gmsh id ptrs for inner circle boundary
            addsuitable_pos_ind_1 = np.where(surf_x_i == mat_connect_additional["pseudo_circle"][:,0])[0] #0 - x
            addsuitable_pos_ind_2 = np.where(surf_y_i == mat_connect_additional["pseudo_circle"][:,1])[0] #1 - y
            addcommon_flag_ptr_i = np.intersect1d(addsuitable_pos_ind_1,addsuitable_pos_ind_2)
            if addcommon_flag_ptr_i.size > 0:
                for addflag_ind_i in range(len(addcommon_flag_ptr_i)):
                    addflag_i = addcommon_flag_ptr_i[addflag_ind_i]
                    additional_mats["pseudo_circle_gmshptrs_connect"][addflag_i, 0] = ptr_i
            addsuitable_pos_ind_3 = np.where(surf_x_i == mat_connect_additional["pseudo_circle"][:,2])[0] #0 - x
            addsuitable_pos_ind_4 = np.where(surf_y_i == mat_connect_additional["pseudo_circle"][:,3])[0] #1 - y
            addcommon_flag_ptr_i = np.intersect1d(addsuitable_pos_ind_3,addsuitable_pos_ind_4)
            if addcommon_flag_ptr_i.size > 0:
                for addflag_ind_i in range(len(addcommon_flag_ptr_i)):
                    addflag_i = addcommon_flag_ptr_i[addflag_ind_i]
                    additional_mats["pseudo_circle_gmshptrs_connect"][addflag_i, 1] = ptr_i

            # #find gmsh id ptrs for outer fault boundary
            # addsuitable_pos_ind_1 = np.where(surf_x_i == mat_connect_additional["fault_outerbc"][:,0])[0] #0 - x
            # addsuitable_pos_ind_2 = np.where(surf_y_i == mat_connect_additional["fault_outerbc"][:,1])[0] #1 - y
            # addcommon_flag_ptr_i = np.intersect1d(addsuitable_pos_ind_1,addsuitable_pos_ind_2)
            # if addcommon_flag_ptr_i.size > 0:
            #     for addflag_ind_i in range(len(addcommon_flag_ptr_i)):
            #         addflag_i = addcommon_flag_ptr_i[addflag_ind_i]
            #         additional_mats["faults_outerbc_gmshptrs_connect"][addflag_i, 0] = ptr_i
            # addsuitable_pos_ind_3 = np.where(surf_x_i == mat_connect_additional["fault_outerbc"][:,2])[0] #0 - x
            # addsuitable_pos_ind_4 = np.where(surf_y_i == mat_connect_additional["fault_outerbc"][:,3])[0] #1 - y
            # addcommon_flag_ptr_i = np.intersect1d(addsuitable_pos_ind_3,addsuitable_pos_ind_4)
            # if addcommon_flag_ptr_i.size > 0:
            #     for addflag_ind_i in range(len(addcommon_flag_ptr_i)):
            #         addflag_i = addcommon_flag_ptr_i[addflag_ind_i]
            #         additional_mats["faults_outerbc_gmshptrs_connect"][addflag_i, 1] = ptr_i

        
    # print(mat_connect_additional["fault_outerbc"])
    # print(additional_mats["faults_outerbc_gmshptrs_connect"])
    # gmsh.model.occ.synchronize()
    # gmsh.fltk.run()
    # exit()
    # np.savetxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv6_cluster/meshv6/matsave.txt",mat_connect_lineendptrs)
    # exit()
    return mat_connect_corner, mat_connect_surf, mat_connect_lineendptrs, additional_mats

def gmsh_add_elems_2nd(mat_connect_corner, 
                       mat_connect_surf, 
                       mat_connect_lineendptrs,
                       additional_mats,
                       case_flag,
                       lc2):

    """
    [inputs]
    mat_connect_corner : connectivity matrix for corner points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
    mat_connect_surf : connectivity matrix for branch points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
    mat_connect_lineendptrs : connect matrix [ptr_i ptr_j ptr_j ptr_k ptr_k ptr_i] such that i-j j-k k-i
    additional_mats : tuple of storing any additional structure 
    case_flag : case specific

    [defines]
    arr_elem_corner : array stores GMSH element id for edge elements   [elem_k]
    arr_elem_surf : array stores GMSH element id for branch elements [elem_k]
    mat_surface_lineloop : matrix for storing surface connectivity [elem1; elem2; elem3; elem4] : elem1 -> elem2 -> elem3 -> elem4 -> elem1
    arr_elem_trialoop : matrix for storing elem id for triangle elements [elem_i elem2 elem3]
    arr_elem_additional : additional elems type
    mat_surface_additional_lineloop : additional surface loop (other than line)

    case_flag = 2D-Cluster
    1-2, 2-3, 3-4, 4-1
    4 - - -[3] - - 3
    |              |
    |              |
   [4]            [2]
    |              |
    1 - - -[1]- -  2

    case_flag: 2D-Cluster-Well
    1-2, 2-3, 3-4, 4-1
    4 - - -[3] - - 3
    |              |
    |              |
    |              |
   [4]      o     [2]
    |              |
    |              |
    |              |
    1 - - [1]- - - 2

    well:
      - - -
    |       |
   [5] [6] [7]
    |       |
      - - - 

    One surface are needed to be created:
    Surface 1: [1] [2] [3] [4]
    Two Loops are needed to be recorded
    """

    #initialize
    arr_elem_corner = np.zeros((np.shape(mat_connect_corner)[0],1),dtype=int)
    arr_elem_surf = np.zeros((np.shape(mat_connect_surf)[0],1),dtype=int)
    arr_elem_trialoop = np.zeros((np.shape(mat_connect_lineendptrs)[0],3),dtype=int)

    #additional surface loop (other than line)
    mat_surface_additional_lineloop = {}
    
    #additional elems type
    arr_elem_additional = {}

    arr_elem_additional["pseudo_circle_gmshlines_connect"] = np.zeros((np.shape(additional_mats["pseudo_circle_gmshptrs_connect"])[0],1))

    #add edge elements 
    for corner_elem_ind in range(np.shape(mat_connect_corner)[0]):
        ptr_i = mat_connect_corner[corner_elem_ind,0]
        ptr_j = mat_connect_corner[corner_elem_ind,1]
        elem_k = gmsh.model.occ.addLine(ptr_i,ptr_j)
        arr_elem_corner[corner_elem_ind] = elem_k

    #add outer square surfacce
    mat_surface_lineloop = np.zeros((4,1),dtype=int)

    ## first surface
    mat_surface_lineloop[0,0] = arr_elem_corner[0];
    mat_surface_lineloop[1,0] = arr_elem_corner[1];
    mat_surface_lineloop[2,0] = arr_elem_corner[2];
    mat_surface_lineloop[3,0] = arr_elem_corner[3];

    loop1 = gmsh.model.occ.addCurveLoop([mat_surface_lineloop[0,0],
                                         mat_surface_lineloop[1,0],       
                                         mat_surface_lineloop[2,0],
                                         mat_surface_lineloop[3,0]])
    
    surf1 = gmsh.model.occ.addPlaneSurface([loop1],tag=10)

    print("finish adding " + str(np.shape(mat_connect_corner)[0]) + " edge elements & surface")

    #-----------------------------------------------------------------------------------------------------------------------------------------#

    track_arc_related_inds = []

    #if it is a circle arc element, do not use addLine, use addCircleArc
    ptr_center = gmsh.model.occ.addPoint(0.0,0.0,0.0,lc2) #always be the center point (0,0)

    #add elem
    for tria_elem_ind in range(np.shape(mat_connect_lineendptrs)[0]):

        #need to find duplicates in this part ......

        flag_arc_related = False

        tria_elem_data = mat_connect_lineendptrs[tria_elem_ind,:]

        if tria_elem_data[0] in additional_mats["pseudo_circle_gmshptrs_connect"] and tria_elem_data[1] in additional_mats["pseudo_circle_gmshptrs_connect"]:
            flag_arc_related = True

        if tria_elem_data[2] in additional_mats["pseudo_circle_gmshptrs_connect"] and tria_elem_data[3] in additional_mats["pseudo_circle_gmshptrs_connect"]:
            flag_arc_related = True

        if tria_elem_data[4] in additional_mats["pseudo_circle_gmshptrs_connect"] and tria_elem_data[5] in additional_mats["pseudo_circle_gmshptrs_connect"]:
            flag_arc_related = True
        
        if flag_arc_related == False:

            elem_k1 = gmsh.model.occ.addLine(tria_elem_data[0],tria_elem_data[1])
            elem_k2 = gmsh.model.occ.addLine(tria_elem_data[2],tria_elem_data[3])
            elem_k3 = gmsh.model.occ.addLine(tria_elem_data[4],tria_elem_data[5])

            arr_elem_trialoop[tria_elem_ind, 0] = elem_k1
            arr_elem_trialoop[tria_elem_ind, 1] = elem_k2
            arr_elem_trialoop[tria_elem_ind, 2] = elem_k3
        
    #arc list
    arc_looplist = []

    #add arc
    for arc_i in range(np.shape(additional_mats["pseudo_circle_gmshptrs_connect"])[0]):

        #get ptr_index
        ptr_i = additional_mats["pseudo_circle_gmshptrs_connect"][arc_i, 0]
        ptr_j = additional_mats["pseudo_circle_gmshptrs_connect"][arc_i, 1]

        #add arc
        arc_k1 = gmsh.model.occ.addCircleArc(int(ptr_i),ptr_center,int(ptr_j))

        #save
        arc_looplist.append(arc_k1)
    
    #inner pseudo circle loop
    loop2 = gmsh.model.occ.addCurveLoop(arc_looplist) 

    surf2 = gmsh.model.occ.addPlaneSurface([loop2],tag=12)

    #cut
    diff1 = gmsh.model.occ.cut([(2,surf1)], [(2, surf2)],removeObject=True)

    # gmsh.model.occ.removeAllDuplicates()
    # gmsh.model.occ.synchronize()
    # gmsh.fltk.run()
    # exit()

    #------------------------------------------------------------------------------------------------------------#
    loop0 = gmsh.model.occ.addCurveLoop([arr_elem_trialoop[0,0],
                                         arr_elem_trialoop[0,1],       
                                         arr_elem_trialoop[0,2]])

    surf0 = gmsh.model.occ.addPlaneSurface([loop0])

    #frag
    fragi = gmsh.model.occ.fragment([(2,surf1)],[(2,surf0)])

    #interior
    num_of_surfaces_nw = np.shape(arr_elem_trialoop)[0]
    for surface_ind in range(1,num_of_surfaces_nw): #43

        if np.sum(arr_elem_trialoop[surface_ind,:]) == 0:
            continue
        else:
            loopi = gmsh.model.occ.addCurveLoop([arr_elem_trialoop[surface_ind,0],
                                                arr_elem_trialoop[surface_ind,1],       
                                                arr_elem_trialoop[surface_ind,2]])

            surfi = gmsh.model.occ.addPlaneSurface([loopi])

            #frag
            fragi = gmsh.model.occ.fragment(fragi[1][0],[(2,surfi)])
    #------------------------------------------------------------------------------------------------------------#

    # entitydim1 = gmsh.model.occ.get_entities(1)
    # print(entitydim1)
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    gmsh.fltk.run()
    # exit()

    return mat_surface_lineloop, arr_elem_corner, arr_elem_surf, arr_elem_trialoop, mat_surface_additional_lineloop, arr_elem_additional

def find_row_index(array, target_row):
    for index in range(np.shape(array)[0]):
        if array[index,:].all() == target_row.all():
            return index
    return -1  # Return -1 if the target row is not found

def gmsh_add_surfaces_2nd(mat_surface_lineloop, 
                          arr_elem_trialoop,
                          mat_surface_additional_lineloop,
                          case_flag):

    """
    [inputs]
    mat_surface_lineloop : matrix for storing surface connectivity [elem1 elem2 elem3 elem4] : elem1 -> elem2 -> elem3 -> elem4 -> elem1
    mat_connect_lineloop : matrix for storing surface network
    case_flag : "2D-Cluster",  "2D-Cluster-Well"
    [defines]
    loopi : line loop i
    surfi : surface i

    case_flag: 2D-Cluster, 2D-Cluster-Well
    1-2, 2-3, 3-4, 4-1
    4 - - -[3] - - 3
    |              |
    |              |
    |              |
   [4]            [2]
    |              |
    |              |
    |              |
    1 - - [1]- - - 2

    physical curves:
    top: 3
    bottom: 1
    left: 4
    right: 2

    additional:
    case_flag: 2D-Cluster-Well

     well:
      - - -
    |       |
   [5] [6] [7]
    |       |
      - - - 

    [returns]
    arr_surface :  array stores GMSH id of surface
    """

    #get number of surfaces
    num_of_surfaces = np.shape(mat_surface_lineloop)[1]

    #initialize arr_surface
    arr_surface = np.zeros((num_of_surfaces,1),dtype=int)

    #additional surface
    arr_additional_surface = {}

    #add surfaces
    for surface_ind in range(num_of_surfaces):
        loopi = gmsh.model.occ.addCurveLoop([mat_surface_lineloop[0,surface_ind],
                                             mat_surface_lineloop[1,surface_ind],       
                                             mat_surface_lineloop[2,surface_ind],
                                             mat_surface_lineloop[3,surface_ind]])
        if case_flag == "2D-Cluster":
            surfi = gmsh.model.occ.addPlaneSurface([loopi])
            arr_surface[surface_ind] = surfi
    
    #add circle surfaces
    if case_flag == "2D-Cluster-Well":
        
        arr_additional_surface["circle"] = np.zeros((np.shape(mat_surface_additional_lineloop["circle"])[1],1),dtype=int)
        
        # line loops #
        #inner circle
        # circle_loopi = gmsh.model.occ.addCurveLoop()

        #store inner boundary (circle-pseudo circle)
        pseudocircle_list = []
        for elem_ind in range(np.shape(mat_surface_additional_lineloop["pseudocircle"])[0]):
            pseudocircle_list.append(mat_surface_additional_lineloop["pseudocircle"][elem_ind,0])
        # print(pseudocircle_list)
        
        #inner pseudo circle loop
        pseudo_circle_loopi = gmsh.model.occ.addCurveLoop(pseudocircle_list) 

        # # surface #   
        # circle_surface = gmsh.model.occ.addPlaneSurface([loopi,circle_loopi],101) #one surface, one circle, hardcode
        # arr_additional_surface["circle"][0] = circle_surface

        pseudo_circle_surface = gmsh.model.occ.addPlaneSurface([pseudo_circle_loopi]) #one surface, one circle, hardcode

        # surface #   
        circle_surface = gmsh.model.occ.addPlaneSurface([loopi]) #one surface, one circle, hardcode
        arr_additional_surface["circle"][0] = circle_surface

    #store network surface
    #get number of surfaces
    num_of_surfaces_nw = np.shape(arr_elem_trialoop)[0]

    #initialize arr_surface
    arr_surface_nw = np.zeros((num_of_surfaces_nw,1),dtype=int)
    #add surfaces
    for surface_ind in range(num_of_surfaces_nw):
        nwloopi = gmsh.model.occ.addCurveLoop([arr_elem_trialoop[surface_ind,0],
                                            arr_elem_trialoop[surface_ind,1],       
                                            arr_elem_trialoop[surface_ind,2]])
        surfi = gmsh.model.occ.addPlaneSurface([nwloopi])
        arr_surface_nw[surface_ind] = surfi
    
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    
    return arr_surface, arr_additional_surface, arr_surface_nw

def gmsh_add_surfaces_occ(mat_surface_lineloop, 
                          arr_elem_trialoop,
                          mat_surface_additional_lineloop,
                          case_flag):
    
    entitydim1 = gmsh.model.occ.get_entities(1)
    # print(entitydim1)
    
    #outer square
    loop1 = gmsh.model.occ.addCurveLoop([mat_surface_lineloop[0,0],
                                         mat_surface_lineloop[1,0],       
                                         mat_surface_lineloop[2,0],
                                         mat_surface_lineloop[3,0]])
    surf1 = gmsh.model.occ.addPlaneSurface([loop1],tag=10)
    
    #inner circle
    # loop2 = gmsh.model.occ.addCurveLoop([mat_surface_additional_lineloop["circle"][0,0],
    #                                      mat_surface_additional_lineloop["circle"][1,0]])
    
    # surf2 = gmsh.model.occ.addPlaneSurface([loop2],tag=11)
    

    #circle channel
    pseudocircle_list = []
    for elem_ind in range(np.shape(mat_surface_additional_lineloop["pseudocircle"])[0]):
        pseudocircle_list.append(mat_surface_additional_lineloop["pseudocircle"][elem_ind,0])

    #inner pseudo circle loop
    loop3 = gmsh.model.occ.addCurveLoop(pseudocircle_list) 

    surf3 = gmsh.model.occ.addPlaneSurface([loop3],tag=12)

    #cut
    diff1 = gmsh.model.occ.cut([(2,surf1)], [(2, surf3)],removeObject=True)

    # #frag
    # # frag1 = gmsh.model.occ.fragment(diff1[1][0], [(2,surf3)])

    entitydim1 = gmsh.model.occ.get_entities(1)
    # print(entitydim1)
    # print(np.shape(entitydim1))
    # exit()

    loop0 = gmsh.model.occ.addCurveLoop([arr_elem_trialoop[0,0],
                                         arr_elem_trialoop[0,1],       
                                         arr_elem_trialoop[0,2]])

    surf0 = gmsh.model.occ.addPlaneSurface([loop0])

    #frag
    fragi = gmsh.model.occ.fragment([(2,surf1)],[(2,surf0)])

    #interior
    num_of_surfaces_nw = np.shape(arr_elem_trialoop)[0]
    for surface_ind in range(1,num_of_surfaces_nw): #43

        # print(surface_ind,num_of_surfaces_nw)
        
        loopi = gmsh.model.occ.addCurveLoop([arr_elem_trialoop[surface_ind,0],
                                             arr_elem_trialoop[surface_ind,1],       
                                             arr_elem_trialoop[surface_ind,2]])

        surfi = gmsh.model.occ.addPlaneSurface([loopi])

        #frag
        fragi = gmsh.model.occ.fragment(fragi[1][0],[(2,surfi)])

    # # #cut
    # diff1 = gmsh.model.occ.cut(fragi[1][0], [(2, surf2)],removeObject=True)

    # #frag
    # frag1 = gmsh.model.occ.fragment(diff1[1][0], [(2,surf3)])

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    entitydim1 = gmsh.model.occ.get_entities(1)
    # print(entitydim1)
    print(np.shape(entitydim1))
    # exit()

    return None


def gmsh_add_physical_entities_2nd(arr_elem_corner, 
                                   arr_elem_surf, 
                                   arr_surface,
                                   arr_elem_additional,
                                   arr_additional_surface,
                                   case_flag):

    """
    [inputs]
    arr_elem_corner : array stores GMSH element id for edge elements   [elem_k]
    arr_elem_surf : array stores GMSH element id for network elements [elem_k]
    arr_surface : array stores GMSH id of surface
    arr_elem_additional : array stores GMSH element id for edge elements   [elem_k]
    arr_additional_surface : array stores GMSH id of non-line-loop surface

    
    case_flag = 2D-Cluster
    4 - - -[3] - - 3
    |              |
    |              |
   [4]             [2]
    |              |
    1 - - -[1]- -  2

    Two surfaces are needed to be created:
    Surface 1: [1] [2] [3] [4]

    physical curves:
    top: 3
    bottom: 1
    left: 4
    right: 2

    additional:
    case_flag: 2D-Cluster-Well

     well:
      - - -
    |       |
   [5] [6] [7]
    |       |
      - - - 

    [returns]
    mat_egelem_phyname: record edge element's physical name [elem_gmsh_id, elem_phy_name]
    mat_nwelem_phyname: record network element's physical name [elem_gmsh_id, elem_phy_name]
    mat_surf_phyname: record surface's physical name [surf_gmsh_id, surf_phy_name] 
    """

    #initialize mat_elem_phyname, mat_surf_phyname
    mat_egelem_phyname = {}
    mat_nwelem_phyname = {}
    mat_surf_phyname = {}

    mat_egelem_phyname_additional = {} #note the key/value pair is different from above
    mat_surface_phyname_additional = {} #note the key/value pair is different from above

    #add boundary physical group
    if case_flag == "2D-Cluster" or "2D-Cluster-Well":
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[2,0]],name="top")
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[0,0]],name="bottom")
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[3,0]],name="left")
        gmsh.model.addPhysicalGroup(1,[arr_elem_corner[1,0]],name="right")

        #save in mat_egelem_phyname
        mat_egelem_phyname[arr_elem_corner[2,0]] = "top"
        mat_egelem_phyname[arr_elem_corner[0,0]] = "bottom"
        mat_egelem_phyname[arr_elem_corner[3,0]] = "left"
        mat_egelem_phyname[arr_elem_corner[1,0]] = "right"
    
    #add circlearc physical group
    if case_flag == "2D-Cluster-Well":
        gmsh.model.addPhysicalGroup(1,[arr_elem_additional["circle"][0,0]],name="well")
        gmsh.model.addPhysicalGroup(1,[arr_elem_additional["circle"][1,0]],name="well")
        mat_egelem_phyname_additional["well"] = [arr_elem_additional["circle"][0,0],arr_elem_additional["circle"][1,0]]
    
    #bf
    # for nw_id in range(np.shape(arr_elem_surf)[0]):
    #     gmsh.model.addPhysicalGroup(1,[arr_elem_surf[nw_id,0]],name="segment" + str(arr_elem_surf[nw_id,0]))
    #     #save in mat_bfelem_phyname
    #     mat_nwelem_phyname[arr_elem_surf[nw_id,0]] = "segment" + str(arr_elem_surf[nw_id,0])

    ##->may need to implement physical segment identification process (remove duplicates)

    # add surface physical group
    if case_flag == "2D-Cluster":
        for surface_id in range(np.shape(arr_surface)[0]):
            gmsh.model.addPhysicalGroup(2,[arr_surface[surface_id,0]],name="surf" + str(arr_surface[surface_id,0]))
            #sabe in mat_surf_phyname
            mat_surf_phyname[arr_surface[surface_id,0]] = "surf" + str(arr_surface[surface_id,0])
    if case_flag == "2D-Cluster-Well":
        for additional_surface_id in range(np.shape(arr_additional_surface["circle"])[0]):
            gmsh.model.addPhysicalGroup(2,[arr_additional_surface["circle"][additional_surface_id,0]],name="surfwell" + str(arr_additional_surface["circle"][additional_surface_id,0]))
            #sabe in mat_surf_phyname
            mat_surface_phyname_additional["surfwell"] = "surfwell" + str(arr_additional_surface["circle"][additional_surface_id,0])

    #synchronize
    gmsh.model.occ.synchronize()

    return mat_egelem_phyname, mat_nwelem_phyname, mat_surf_phyname, mat_egelem_phyname_additional, mat_surface_phyname_additional


#generate SubDomian ID
def meshio_generate_subdomainID(file_path2,
                                txt_elems_path,
                                txt_ids_path,
                                txt_surround_path,
                                mat_surf_triacoordall):

    """
    [inputs]
    file_path2 : file path for saving second meshing file
    txt_elems_path,txt_ids_path,txt_surround_path,
    mat_surf_tiracoordall : save coordinates for three points of a tria elem
    
    Note : [ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y, ptr_k_x, ptr_k_y]
    """

    #read file
    m = meshio.read(file_path2)

    #Get TRIA3 Connectivity
    tria_elem_connect = m.cells_dict['triangle']

    #TRIA3 Elem Num
    num_elem = np.shape(tria_elem_connect)[0]

    print("total elems:",num_elem)

    #Define arr object for storing centroid coordinates for each coarse tria element
    arr_centriod_tria = np.zeros((np.shape(mat_surf_triacoordall)[0],2))
    arr_centriod_tria[:,0] = (mat_surf_triacoordall[:,0] + mat_surf_triacoordall[:,2] + mat_surf_triacoordall[:,4])/3
    arr_centriod_tria[:,1] = (mat_surf_triacoordall[:,1] + mat_surf_triacoordall[:,3] + mat_surf_triacoordall[:,5])/3

    #Define dict object for storing subdomain elements
    arr_index_elem = {}

    #Loop over elem
    for elem_ind in range(num_elem):

        print("elem_ind: ", elem_ind)
        print("Progress: ", elem_ind/num_elem*100,"%")
        
        #get connectivity of current element
        elem_connect_i = tria_elem_connect[elem_ind,:]

        #get end points coordinates for current element
        coord_data_x = m.points[elem_connect_i][:,0]
        coord_data_y = m.points[elem_connect_i][:,1]

        #get centroid point coordinate for current element
        coord_centriod_x = np.sum(coord_data_x) / 3
        coord_centriod_y = np.sum(coord_data_y) / 3

        #create arr contains the centriod ptrs
        arr_centriod_duplicate = np.zeros((np.shape(mat_surf_triacoordall)[0],2))
        arr_centriod_duplicate[:,0] = coord_centriod_x
        arr_centriod_duplicate[:,1] = coord_centriod_y

        #find the distance 
        dist_mat = distance.cdist(arr_centriod_tria,arr_centriod_duplicate)
        
        #sort the column from min to max, get corresponding index
        dist_mat = dist_mat[:,0]
        ind_sort = np.argsort(dist_mat)

        for ind in range(len(ind_sort)):

            #get current element
            tria_elem_i = mat_surf_triacoordall[ind_sort[ind]]
            
            #get coordinate in x and y
            x1 = tria_elem_i[0]; y1 = tria_elem_i[1];
            x2 = tria_elem_i[2]; y2 = tria_elem_i[3];
            x3 = tria_elem_i[4]; y3 = tria_elem_i[5];

            #check whether point lies in the triangle
            bool_isinside = func_isInside(x1,y1,x2,y2,x3,y3,coord_centriod_x,coord_centriod_y)

            #save pair: tria_elem_ind : elem_ind
            if bool_isinside:
                
                #if no elems defined before, create a list contains the current elem
                if arr_index_elem.get(ind_sort[ind]) == None:
                    arr_index_elem[ind_sort[ind]] = [elem_ind]
                else: #if already defined, append the current value to existing list
                    arr_index_elem[ind_sort[ind]].append(elem_ind)
                
                #break for loop
                break
        
    #loop over elems in dict
    #define list_subdomain_elem
    list_subdomain_elem = []
    
    #define list_subdomain_ids
    list_subdomain_ids = []

    list_subdomain_surrounding = []
    
    #loop over keys
    list_keys = list(arr_index_elem.keys())
    for key_i in list_keys:
        
        #get list of values
        list_vals_i = arr_index_elem[key_i]

        #create subdomain ids
        list_ids_i = [key_i + 1] * len(list_vals_i)

        #save 
        list_subdomain_elem.extend(list_vals_i)
        list_subdomain_ids.extend(list_ids_i)
        list_subdomain_surrounding.append(key_i+1)
    
    #save txt file
    np.savetxt(txt_elems_path,list_subdomain_elem,fmt='%i',newline=" ")
    np.savetxt(txt_ids_path,list_subdomain_ids,fmt='%i',newline=" ")
    np.savetxt(txt_surround_path,list_subdomain_surrounding,fmt='%i',newline=" ")

def meshio_find_physical_groups(file_path2,
                                txt_global_path,
                                array_coord,
                                radius_2nd,
                                lc2):
    
    #read file
    m = meshio.read(file_path2)

    #Get TRIA3 Connectivity
    tria_elem_connect = m.cells_dict['triangle']

    #Get coords
    bound_min_x = array_coord[0,0] #x1
    bound_max_x = array_coord[1,0] #x2
    bound_min_y = array_coord[0,1] #y1
    bound_max_y = array_coord[2,1] #y3

    borehole_center_x = array_coord[5,0]
    borehole_center_y = array_coord[5,1]

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

def func_isInside(x1,y1,x2,y2,x3,y3,coord_centriod_x,coord_centriod_y):

    x_A = x1; y_A = y1
    x_B = x2; y_B = y2
    x_C = x3; y_C = y3
    x_P = coord_centriod_x
    y_P = coord_centriod_y

    new_x_B = x_B - x_A
    new_y_B = y_B - y_A
    new_x_C = x_C - x_A 
    new_y_C = y_C - y_A
    new_x_P = x_P - x_A
    new_y_p = y_P - y_A

    w_A_upper = new_x_P * ( new_y_B - new_y_C ) + new_y_p * ( new_x_C - new_x_B ) + new_x_B * new_y_C - new_x_C * new_y_B
    w_B_upper = new_x_P * new_y_C - new_y_p * new_x_C
    w_C_upper = new_y_p * new_x_B - new_x_P * new_y_B

    d = new_x_B * new_y_C - new_x_C * new_y_B

    if d  > 1e-13:
        w_A = w_A_upper / d
        w_B = w_B_upper / d
        w_C = w_C_upper / d
    else:
        w_A = 0
        w_B = 0
        w_C = 0
    
    if ( w_A >= 0.0 and w_A <= 1.0 and w_B >= 0.0 and w_B <= 1.0 and w_C >= 0.0 and w_C <= 1.0 ):
        return True
    else:
        return False

# computes the alpha-shape (concave hull) and keeps only the outer boundary.
# reference: https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points 
def helper_alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

if __name__ == "__main__":

    from matplotlib.pyplot import *

    #testcode : https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    # Constructing the input point data
    # np.random.seed(0)
    # x = 3.0 * np.random.rand(2000)
    # y = 2.0 * np.random.rand(2000) - 1.0
    # inside = ((x ** 2 + y ** 2 > 1.0) & ((x - 3) ** 2 + y ** 2 > 1.0))
    # points = np.vstack([x[inside], y[inside]]).T

    #own data
    faultpoints = np.loadtxt("/Users/andyz/projects/febe_moose/problems_febe/test_CZM_multifaultsv6_cluster/meshv6/faultpoints.txt")
    points = np.array(faultpoints).reshape([int(len(faultpoints)/2),2])
    print(points)

    # Computing the alpha shape
    edges = helper_alpha_shape(points, alpha=800, only_outer=True)

    print(edges)
    print(len(edges))

    #get unique points
    ptrs_outerbc_list = []
    for i, j in edges:
        ptrs_outerbc_list.append(i)
        ptrs_outerbc_list.append(j)
    arr_ptrs_outerbc = np.unique(ptrs_outerbc_list)

    print(arr_ptrs_outerbc)

    # get points coordinates
    points_outerbc = points[arr_ptrs_outerbc,:]

    # zero as initial points
    arr_origin_x = np.zeros(np.shape(points_outerbc)[0])
    arr_origin_y = np.zeros(np.shape(points_outerbc)[0])
    #determine sequence #hardcode origin (0,0)
    angles = helper_angle_between_points(points_outerbc[:,0], points_outerbc[:,1], arr_origin_x, arr_origin_y)
    #sort 
    sorted_ind = np.argsort(angles)
    outerbc_x_sorted = points_outerbc[:,0][sorted_ind]
    outerbc_y_sorted = points_outerbc[:,1][sorted_ind]
    #add initial one at end
    outerbc_x_sorted = np.append(outerbc_x_sorted,outerbc_x_sorted[0])
    outerbc_y_sorted = np.append(outerbc_y_sorted,outerbc_y_sorted[0])

    # Generate line segement matrix [ptri_x ptri_y ptrj_x ptrj_y] line segment i -> j
    mat_connect_additional = {}
    mat_connect_additional["fault_outerbc"] = np.zeros((len(edges),4))
    for line_i in range(np.shape(points_outerbc)[0]):
        mat_connect_additional["fault_outerbc"][line_i,0] = outerbc_x_sorted[line_i]
        mat_connect_additional["fault_outerbc"][line_i,1] = outerbc_y_sorted[line_i]
        mat_connect_additional["fault_outerbc"][line_i,2] = outerbc_x_sorted[line_i+1]
        mat_connect_additional["fault_outerbc"][line_i,3] = outerbc_y_sorted[line_i+1]
    
    print(mat_connect_additional["fault_outerbc"])

    # Plotting the output
    figure()
    axis('equal')
    plot(points[:, 0], points[:, 1], '.')
    for i, j in edges:
        plot(points[[i, j], 0], points[[i, j], 1])
    show()