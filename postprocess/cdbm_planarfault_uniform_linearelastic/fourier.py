#import modulus
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import scipy

#-----------------------------------------------------functions-----------------------------------------------------#
def DecodeName(exodus_file_path, decodeflag,save_folder_output_file_path):

    """
    Get Element Var Name 
    """

    #Read File
    nc = netCDF4.Dataset(exodus_file_path)

    #Get numpy.bytes name array
    if decodeflag == "name_elem_var":
        arr_name_elem_var = nc.variables['name_elem_var']
    elif decodeflag == "name_nod_var":
        arr_name_elem_var = nc.variables['name_nod_var']
    elif decodeflag == "eb_names":
        arr_name_elem_var = nc.variables['eb_names']

    #Get number of name
    num_names = np.shape(arr_name_elem_var)[0]

    #initialize list of names
    list_name_elem_var = []

    #loop over each name
    for name_ind in range(num_names):

        name_i_decode = ''

        name_i = arr_name_elem_var[name_ind]

        # print(name_i)

        for ind in range(len(name_i)):
            
            if not np.ma.is_masked(name_i[ind]):
                
                byte_i_decode = name_i[ind].decode('UTF-8')

                name_i_decode += byte_i_decode

        list_name_elem_var.append(name_i_decode)
    
    #save
    if decodeflag == "name_elem_var":
        np.savetxt(save_folder_output_file_path + "/list_name_elem_var.txt",list_name_elem_var,fmt='%s',newline=" ")
    elif decodeflag == "name_nod_var":
        np.savetxt(save_folder_output_file_path + "/list_name_nod_var.txt",list_name_elem_var,fmt='%s',newline=" ")
    elif decodeflag == "eb_names":
        np.savetxt(save_folder_output_file_path + "/list_eb_names.txt",list_name_elem_var,fmt='%s',newline=" ")
#-------------------------------------------------------------------------------------------------------------------#

#exodus file
exodus_file_path = "/Volumes/One Touch/Research/DamageBreakage/planarfaulttests/accel_fourier_transform/test_planarfault_main_out_linear_elastic.e"

#save path
save_path = "./outputs/fourier"

#velocity (temp)
plot_var_name_0 = "accel_slipweakening_x"
plot_var_name_1 = "accel_slipweakening_y"

#point coord (variable)
ptr_x = 6019.45
ptr_y = -701.224

#read
nc = netCDF4.Dataset(exodus_file_path)

#decode element names
DecodeName(exodus_file_path,"name_nod_var",save_path)

#Initialize dict storing elem var name / index pair
dict_evn_index = {}

#Read Elem Var Name
evn = np.loadtxt(save_path + "/list_name_nod_var.txt",dtype=str)

#fill in dict_evn_index
for evn_ind in range(len(evn)):
    evn_i = evn[evn_ind]
    # start from 1
    dict_evn_index[evn_i] = evn_ind + 1

var_ind_0 = dict_evn_index[plot_var_name_0]
var_ind_1 = dict_evn_index[plot_var_name_1]

arr_accel_x = nc.variables["vals_nod_var"+str(var_ind_0)]
arr_accel_y = nc.variables["vals_nod_var"+str(var_ind_1)]

arr_xcoord = nc.variables['coordx'][:]
arr_ycoord = nc.variables['coordy'][:]

ptr_index_x = np.where(abs(arr_xcoord - ptr_x) < 1e-2)[0]
ptr_index_y = np.where(abs(arr_ycoord - ptr_y) < 1e-2)[0]
ptr_index = np.intersect1d(ptr_index_x,ptr_index_y)[0]

print("The point coordinate: ", arr_xcoord[ptr_index],arr_ycoord[ptr_index])

#get data for this point
arr_accel_x_ptr = arr_accel_x[:, ptr_index]
arr_accel_y_ptr = arr_accel_y[:, ptr_index]

#get mag
arr_accel_mag_ptr = np.sqrt(arr_accel_x_ptr ** 2 + arr_accel_y_ptr ** 2)

#time
arr_time = nc.variables['time_whole'][:]

accelfft = np.fft.fft(arr_accel_mag_ptr[:len(arr_time)])

tstep = 0.05             #sampling time interval
Fs = 1 / tstep           #sampling freq
N = len(accelfft)        #number of samples
n = np.arange(N)         #sample arr
T = N/Fs                 #total time
freq = n/T               #freq arr

accelfft_normalized = np.abs(accelfft) / N

#
freq_plot = freq[0:int(N/2+1)]
accelfft_plot = 2 * accelfft_normalized[0:int(N/2+1)]
accelfft_plot[0] = accelfft_plot[0] / 2

#normalize by 0 freq value
normalized_ratio = accelfft_plot[1]

accelfft_plot = accelfft_plot / normalized_ratio

plt.figure()
plt.loglog(freq_plot,accelfft_plot,'ro-',label='intact')
plt.xlabel("Freq (Hz)")
plt.ylabel("Accel")
# plt.legend()
# plt.show()

# exit()
#-------------------------------------------------------------------------------------------------------------------#

#exodus file
exodus_file_path = "/Volumes/One Touch/Research/DamageBreakage/planarfaulttests/accel_fourier_transform/test_planarfault_main_out_damage.e"

#save path
save_path = "./outputs/fourier"

#velocity (temp)
plot_var_name_0 = "accel_slipweakening_x"
plot_var_name_1 = "accel_slipweakening_y"

# #point coord (variable)
# ptr_x = -3432.07
# ptr_y = 353.977

#read
nc = netCDF4.Dataset(exodus_file_path)

#decode element names
DecodeName(exodus_file_path,"name_nod_var",save_path)

#Initialize dict storing elem var name / index pair
dict_evn_index = {}

#Read Elem Var Name
evn = np.loadtxt(save_path + "/list_name_nod_var.txt",dtype=str)

#fill in dict_evn_index
for evn_ind in range(len(evn)):
    evn_i = evn[evn_ind]
    # start from 1
    dict_evn_index[evn_i] = evn_ind + 1

var_ind_0 = dict_evn_index[plot_var_name_0]
var_ind_1 = dict_evn_index[plot_var_name_1]

arr_accel_x = nc.variables["vals_nod_var"+str(var_ind_0)]
arr_accel_y = nc.variables["vals_nod_var"+str(var_ind_1)]

arr_xcoord = nc.variables['coordx'][:]
arr_ycoord = nc.variables['coordy'][:]

ptr_index_x = np.where(abs(arr_xcoord - ptr_x) < 1e-2)[0]
ptr_index_y = np.where(abs(arr_ycoord - ptr_y) < 1e-2)[0]
ptr_index = np.intersect1d(ptr_index_x,ptr_index_y)[0]

print("The point coordinate: ", arr_xcoord[ptr_index],arr_ycoord[ptr_index])

#get data for this point
arr_accel_x_ptr = arr_accel_x[:, ptr_index]
arr_accel_y_ptr = arr_accel_y[:, ptr_index]

#get mag
arr_accel_mag_ptr = np.sqrt(arr_accel_x_ptr ** 2 + arr_accel_y_ptr ** 2)

#time
arr_time = nc.variables['time_whole'][:]

accelfft = np.fft.fft(arr_accel_mag_ptr[:len(arr_time)])

tstep = 0.05             #sampling time interval
Fs = 1 / tstep           #sampling freq
N = len(accelfft)        #number of samples
n = np.arange(N)         #sample arr
T = N/Fs                 #total time
freq = n/T               #freq arr

accelfft_normalized = np.abs(accelfft) / N

#
freq_plot2 = freq[0:int(N/2+1)]
accelfft_plot2 = 2 * accelfft_normalized[0:int(N/2+1)]
accelfft_plot2[0] = accelfft_plot2[0] / 2

#normalize by 0 freq value
normalized_ratio = accelfft_plot2[1]

accelfft_plot2 = accelfft_plot2 / normalized_ratio

# print(accelfft_plot2[1], accelfft_plot[1])

# plt.figure()
plt.loglog(freq_plot2,accelfft_plot2,'bo-',label='damage')
plt.xlabel("Freq (Hz)")
plt.ylabel("Accel")
plt.legend()
plt.show()

# exit()