import netCDF4
import numpy as np

def DecodeName(nc, decodeflag, save_folder_output_file_path):

    """
    Get Nodal Var Name 
    """

    # #Read File
    # nc = netCDF4.Dataset(exodus_file_path)

    #Get numpy.bytes name array
    if decodeflag == "name_elem_var":
        arr_name_nod_var = nc.variables['name_elem_var']

    #Get number of name
    num_names = np.shape(arr_name_nod_var)[0]

    #initialize list of names
    list_name_elem_var = []

    #loop over each name
    for name_ind in range(num_names):

        name_i_decode = ''

        name_i = arr_name_nod_var[name_ind]

        # print(name_i)

        for ind in range(len(name_i)):
            
            if not np.ma.is_masked(name_i[ind]):
                
                byte_i_decode = name_i[ind].decode('UTF-8')

                name_i_decode += byte_i_decode

        list_name_elem_var.append(name_i_decode)
    
    #save
    if decodeflag == "name_elem_var":
        np.savetxt(save_folder_output_file_path + "/list_name_elem_var.txt",list_name_elem_var,fmt='%s',newline=" ")

    return list_name_elem_var


if __name__ == '__main__':

    #file path
    exodus_file_path = "/Users/zhaoc/Downloads/tpv2053D/tpv2053D_out_jun18.e"
    save_folder_output_file_path = "./farms_data"

    #read exodus file
    nc = netCDF4.Dataset(exodus_file_path)

    #decode name
    decodeflag = "name_elem_var"

    ##Decode name
    elemvarnames = DecodeName(nc, decodeflag ,save_folder_output_file_path)

    #blockid (which contains czm data)
    elemblock_ind = 2

    #get connectivity 
    tet_elem_connect = nc.variables['connect' + str(elemblock_ind)][:]-1

    #tet4 Elem Num
    num_elem = np.shape(tet_elem_connect)[0]

    #define list
    #initialize list for saving elem along main fault
    list_elem_coords = []

    #loop over element
    for elem_ind in range(num_elem):
    
        #get connectivity of current element
        elem_connect_i = tet_elem_connect[elem_ind,:]

        #get coordinate for current element
        coord_data_x = nc.variables['coordx'][elem_connect_i]
        coord_data_y = nc.variables['coordy'][elem_connect_i]
        coord_data_z = nc.variables['coordz'][elem_connect_i]

        #get average (central ptr coordinates)
        avg_coordx = np.sum(coord_data_x)/len(coord_data_x)
        avg_coordy = np.sum(coord_data_y)/len(coord_data_y)
        avg_coordz = np.sum(coord_data_z)/len(coord_data_z)

        #save elem info
        list_elem_coords.extend([elem_ind, avg_coordx, avg_coordy, avg_coordz])

    print(nc.variables.keys())