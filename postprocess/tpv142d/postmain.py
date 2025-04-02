import sys
import os

#this path needs to be changed everytime you put the main file in a different folder
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import postprocess.postfunc as postfunc

#parameters
plotvar = ["jump_x","jump_x_rate","traction_x"]
save_file_path = "./results"
dim = 2
node_per_elem = 3 #hex8 #4 #tet4
elemblock_ind = 2 #primary surface associated blocks

#clear folders before postprocess
folders_remove_path = './'
folders_to_remove = ['farms_data_elemental', 'outputs']

#whether to save the folder results
rename_folder = True
rename_name = "./farms_data_elemental"

#multiple nodes
ptr_coords = [[ 4500 , 0, -7500],
              [-7500 , 0, -7500],
              [-12000, 0, -7500],
              [   0  , 0, -7500],              
              [4500  , 0, -7500],
              [7500  , 0, -7500],
              [12000 , 0, -7500]]

#multiple files
decodeflags = ["name_elem_var"]
# file_paths = ["/Users/chunhuizhao/Downloads/tpv2052d/tpv2052D_quad_out.e"]
# file_paths = ["/Users/chunhuizhao/projects/farms_benchmark/examples/benchmark_tpv2052D/czm/tpv2052D_quad_out.e"]
# file_paths = ["/Users/chunhuizhao/projects/farms_benchmark/examples/benchmark_tpv2052D/czm/tpv2052D_tria_out.e"]
file_paths = ["/Users/chunhuizhao/projects/farms_benchmark/examples/benchmark_tpv2052D/czm/tpv2052D_tria_directcentral_out.e"]

if __name__ == '__main__':

    #ensure the save path folder exists
    postfunc.systemops.ensure_folder_exists(save_file_path)

    #remove all files in results folder
    postfunc.systemops.remove_specific_folders(folders_remove_path, folders_to_remove)

    #read exodus file, initialize variables
    ppc = postfunc.postprocessclass(file_path=file_paths[0], 
                                    decodeflags=decodeflags, 
                                    plotvar=plotvar, 
                                    save_file_path=save_file_path,
                                    dim=dim,
                                    node_per_elem=node_per_elem)

    #decode the element & nodal properties
    ppc.decode_name()

    #pre-post get element info
    ppc.pre_post(elemblock_ind = elemblock_ind)
        
    #loop over multiple nodes
    for j in range(len(ptr_coords)):

        #get ptr_coord
        ptr_coord = ptr_coords[j]

        #calculate nodal point quantities
        ppc.post_elemental(ptr_coord=ptr_coord,
                            elemblock_ind = elemblock_ind)
    
    #rename the folder
    if rename_folder:
        postfunc.systemops.rename_folder(save_file_path, rename_name)
