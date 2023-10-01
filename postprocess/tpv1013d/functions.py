import numpy as np

def DecodeNameNodal(nc, decodeflag, save_folder_output_file_path):

    """
    Get Nodal Var Name 
    """

    # #Read File
    # nc = netCDF4.Dataset(exodus_file_path)

    #Get numpy.bytes name array
    if decodeflag == "name_nod_var":
        arr_name_nod_var = nc.variables['name_nod_var']

    #Get number of name
    num_names = np.shape(arr_name_nod_var)[0]

    #initialize list of names
    list_name_nod_var = []

    #loop over each name
    for name_ind in range(num_names):

        name_i_decode = ''

        name_i = arr_name_nod_var[name_ind]

        # print(name_i)

        for ind in range(len(name_i)):
            
            if not np.ma.is_masked(name_i[ind]):
                
                byte_i_decode = name_i[ind].decode('UTF-8')

                name_i_decode += byte_i_decode

        list_name_nod_var.append(name_i_decode)
    
    #save
    if decodeflag == "name_nod_var":
        np.savetxt(save_folder_output_file_path + "/list_name_nod_var.txt",list_name_nod_var,fmt='%s',newline=" ")