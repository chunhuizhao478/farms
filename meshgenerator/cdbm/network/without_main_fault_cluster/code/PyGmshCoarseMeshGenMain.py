"""
Coarser Mesh Generation

Created: Chunhui Zhao, Aug 14, 2023
"""

# import module 
import numpy as np
import gmsh
import math
import matplotlib.pyplot as plt
from scipy.stats import qmc
import meshio
from scipy.spatial import distance

#import function file
import PyGmshBuildExplodeFunc as efunc

# First Step Mesh Size 
lc1 = 2000

# Meshing Coarse Mesh Size
lc2 = 800 #for outer domain
lc3 = 25  #for inner network region

#list_cornerptrs: list of corner ptrs coordinates [x1, y1, x2, y2, x3, y3, x4, y4]
list_cornerptrs = [-16000.0,-16000.0,
                    16000.0,-16000.0,
                    16000.0,16000.0,
                   -16000.0,16000.0]

# Region of Multifaults
dist_verticalbc = 6500
dist_lateralbc = 6500

# Meshing Algorithm
num_alg_meshing = 5

# global path
global_path = "/Users/andyz/projects/farms/meshgenerator/cdbm/network/without_main_fault_cluster"

##file path
#file saving coarse mesh
file_path1 = global_path + "/mshfiles/coarse/multifaults_coarse.msh"
file_path2 = global_path + "/mshfiles/coarse/network_coarse.msh"
txt_elems_path = global_path + "/elemdata/coarse/elem_ind.txt"
txt_ids_path = global_path + "/elemdata/coarse/elem_ids.txt"
txt_surround_path = global_path + "/elemdata/coarse/elem_surroundingblocks.txt"
txt_block_pairs = global_path + "/elemdata/coarse/elem_blockpairs_bruteforce.txt"

csv_file_path = global_path + "/elemdata/coarse/area.csv"

png1_file_path = global_path + "/elemdata/coarse/1stptr.png"
png2_file_path = global_path + "/elemdata/coarse/lineplot.png"

#
txt_mat_surf_ptr_all = global_path + "/elemdata/coarse/mat_surf_ptr_all.txt"
txt_mat_surf_connect_all = global_path + "/elemdata/coarse/mat_surf_connect_all.txt"
txt_mat_surf_triacoordall = global_path + "/elemdata/coarse/mat_surf_triacoordall.txt"
txt_mat_surf_trialinesegmentall = global_path + "/elemdata/coarse/mat_surf_trialinesegmentall.txt"


#
case_flag = "2D-Cluster"
radius = None

######## Generation Begins ########

##------Mesh on Coarse Grid------##
#return array of corner coordinates [coord_x,coord_y](numofnodes,2)
array_coord = efunc.gen_coord_1st(list_cornerptrs=list_cornerptrs,
                                  case_flag=case_flag)

#return connectivity of boundary lines [ptr_i,ptr_j] (numoflines,2)
#return additional_mats : tuple of storing any additional structure 
#additional_mats["circle"] for access
mat_connect_corner, additional_mats = efunc.gmsh_add_ptrs_1st(array_coord=array_coord,
                                                                lc1 = lc1,
                                                                case_flag=case_flag)

#return line loop, elemid of outer boundary
#return line loop, elemid of other non-line boundary
mat_surface_lineloop, arr_elem_corner,mat_surface_additional_lineloop, arr_elem_additional = efunc.gmsh_add_elems_1st(mat_connect_corner=mat_connect_corner,
                                                                                                                        additional_mats=additional_mats,
                                                                                                                        case_flag=case_flag)

#return array stores GMSH id of surface
arr_surface,arr_additional_surface = efunc.gmsh_add_surfaces_1st(mat_surface_lineloop=mat_surface_lineloop,
                                                                 mat_surface_additional_lineloop=mat_surface_additional_lineloop,
                                                                 case_flag=case_flag)
print(arr_additional_surface)

#return edge, surface physical name
mat_egelem_phyname, mat_surf_phyname, mat_egelem_phyname_additional, mat_surface_phyname_additional = efunc.gmsh_add_physical_entities_1st(arr_elem_corner=arr_elem_corner, 
                                                                                                                                            arr_elem_additional=arr_elem_additional,
                                                                                                                                            arr_surface=arr_surface,
                                                                                                                                            arr_additional_surface=arr_additional_surface,
                                                                                                                                            case_flag=case_flag)

print(mat_egelem_phyname, mat_surf_phyname, mat_egelem_phyname_additional, mat_surface_phyname_additional)

#return coarse mesh file
efunc.gmsh_meshing(num_alg_meshing,file_path1) 

#return 
#mat_surf1_connect: Define Larger mat for store fault line connectivity [ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y]
#mat_surf1_triacoordall: Define Larger mat for store triangle three coordinates [ptr_i_x, ptr_i_y, ptr_j_x, ptr_j_y, ptr_k_x, ptr_k_y]
#mat_surf_ptr_all: Define Larger mat for store coordinates [ptr_i_x, ptr_i_y]
mat_surf_ptr_all, mat_surf_connect_all, mat_surf_triacoordall, mat_surf_trialinesegmentall, mat_ptrs_additional, mat_connect_additional = efunc.meshio_save_network_info(file_path1=file_path1,
                                                                                                                                                                        array_coord=array_coord,
                                                                                                                                                                        dist_verticalbc=dist_verticalbc,
                                                                                                                                                                        dist_lateralbc=dist_lateralbc,
                                                                                                                                                                        case_flag=case_flag,
                                                                                                                                                                        surfnames = mat_surf_phyname,
                                                                                                                                                                        surfnames_additional = mat_surface_phyname_additional,
                                                                                                                                                                        radius=radius)

#save coarse mesh info
np.savetxt(txt_mat_surf_ptr_all,mat_surf_ptr_all)
np.savetxt(txt_mat_surf_connect_all,mat_surf_connect_all)
np.savetxt(txt_mat_surf_triacoordall,mat_surf_triacoordall)
np.savetxt(txt_mat_surf_trialinesegmentall,mat_surf_trialinesegmentall)

##------Mesh on Coarse Grid End------##

##------Mesh on Fine Grid------##

array_coord = efunc.helper_modify_array_coord(array_coord=array_coord,
                                              list_circleptrs_2nd=None,
                                              case_flag=case_flag)

#return 
#mat_connect_corner : connectivity matrix for corner points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
#mat_connect_surf : connectivity matrix for branch points [ptr_i, ptr_j] : ptr_i -> ptr_j (arrow indicates the direction)
#additional_mats : tuple of storing any additional structure 
mat_connect_corner, mat_connect_surf, mat_connect_lineendptrs, additional_mats =  efunc.gmsh_add_points_2nd(array_coord=array_coord, 
                                                                                                            mat_surf_ptr=mat_surf_ptr_all,
                                                                                                            mat_surf_connect=mat_surf_connect_all,
                                                                                                            mat_surf_trialinesegmentall=mat_surf_trialinesegmentall,
                                                                                                            mat_connect_additional = mat_connect_additional,
                                                                                                            lc2=lc2,
                                                                                                            lc3=lc3,
                                                                                                            case_flag=case_flag)
#return
#arr_elem_corner : array stores GMSH element id for edge elements   [elem_k]
#arr_elem_surf : array stores GMSH element id for branch elements [elem_k]
#mat_surface_lineloop : matrix for storing surface connectivity [elem1; elem2; elem3; elem4] : elem1 -> elem2 -> elem3 -> elem4 -> elem1
#arr_elem_trialoop : matrix for storing elem id for triangle elements [elem_i elem2 elem3]
#mat_surface_additional_lineloop : additional surface loop (other than line)
#arr_elem_additional : additional elems type
mat_surface_lineloop, arr_elem_corner, arr_elem_surf, arr_elem_trialoop, mat_surface_additional_lineloop, arr_elem_additional = efunc.gmsh_add_elems_2nd(mat_connect_corner=mat_connect_corner, 
                                                                                                                                                         mat_connect_surf=mat_connect_surf, 
                                                                                                                                                         mat_connect_lineendptrs=mat_connect_lineendptrs,
                                                                                                                                                         additional_mats=additional_mats,
                                                                                                                                                         case_flag=case_flag)

#return
#arr_surface :  array stores GMSH id of surface
#arr_additional_surface :  array stores GMSH id of non-planar surface
#arr_surface_nw :  array stores GMSH id of fault-triangle surface
arr_surface, arr_additional_surface, arr_surface_nw = efunc.gmsh_add_surfaces_2nd(mat_surface_lineloop=mat_surface_lineloop, 
                                                                                  arr_elem_trialoop=arr_elem_trialoop,
                                                                                  mat_surface_additional_lineloop=mat_surface_additional_lineloop,
                                                                                  case_flag=case_flag)

# #Test #works for well
# efunc.gmsh_add_surfaces_occ(mat_surface_lineloop=mat_surface_lineloop, 
#                             arr_elem_trialoop=arr_elem_trialoop,
#                             mat_surface_additional_lineloop=mat_surface_additional_lineloop,
#                             case_flag=case_flag)


#return
#arr_elem_corner : array stores GMSH element id for edge elements   [elem_k]
#arr_elem_surf : array stores GMSH element id for network elements [elem_k]
#arr_surface : array stores GMSH id of surface
mat_egelem_phyname, mat_nwelem_phyname, mat_surf_phyname, mat_egelem_phyname_additional, mat_surface_phyname_additional = efunc.gmsh_add_physical_entities_2nd(arr_elem_corner=arr_elem_corner, 
                                                                                                                                                                arr_elem_surf=arr_elem_surf, 
                                                                                                                                                                arr_surface=arr_surface,
                                                                                                                                                                arr_elem_additional=arr_elem_additional,
                                                                                                                                                                arr_additional_surface=arr_additional_surface,
                                                                                                                                                                case_flag=case_flag)

# #add boxes
# #gmsh add bounding box
# gmsh.model.mesh.field.add("Box", 1)
# gmsh.model.mesh.field.setNumber(1, "VIn", lc2/8)
# gmsh.model.mesh.field.setNumber(1, "VOut", lc2)
# gmsh.model.mesh.field.setNumber(1, "XMin", -11000)
# gmsh.model.mesh.field.setNumber(1, "XMax", 11000)
# gmsh.model.mesh.field.setNumber(1, "YMin", -11000)
# gmsh.model.mesh.field.setNumber(1, "YMax", 11000)
# gmsh.model.mesh.field.setNumber(1, "Thickness", 2000)

# gmsh.model.mesh.field.add("Min", 3)
# gmsh.model.mesh.field.setNumbers(3, "FieldsList", [1])
# gmsh.model.mesh.field.setAsBackgroundMesh(3)

#return final mesh file 
efunc.gmsh_meshing(num_alg_meshing, file_path2)

#return saved subdomainIDs & elementIDs
efunc.meshio_generate_subdomainID(file_path2=file_path2,
                                  txt_elems_path=txt_elems_path,
                                  txt_ids_path=txt_ids_path,
                                  txt_surround_path=txt_surround_path,
                                  txt_block_pairs = txt_block_pairs,
                                  mat_surf_triacoordall=mat_surf_triacoordall)