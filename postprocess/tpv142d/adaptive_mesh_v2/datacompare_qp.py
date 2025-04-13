import numpy as np
import matplotlib.pyplot as plt
import os

#benchmark code
benchmark_code = "TPV14-2D"

#read benchmark data
benchmark_label = "benchmark-DG-200m"

#read farms data
farms_label = "farms-100m-qp"

#time farms 0.1s interval
time = np.linspace(0.050,12.0,240)

#check dir
def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory '{directory}' created.")
    else:
        print(f"Directory '{directory}' already exists.")

def find_and_read_file(filename, search_directory):
    # Walk through the directory
    for root, dirs, files in os.walk(search_directory):
        if filename in files:
            file_path = os.path.join(root, filename)
            print(f"File found: {file_path}")
            with open(file_path, 'r') as file:
                contents = file.read()
            return contents
    
    # If the file is not found
    print(f"File '{filename}' not found in directory '{search_directory}'.")
    return None

#define plot function
def plotfigure(benchmark_data, time, farms_sliprate_data, farms_slip_data, farms_traction_data, farms_traction_y, benchmark_label, farms_label, strike_value, dip_value, benchmark_code):

    ## check if dir exists, create one if not
    saved_path = "./outputs_qp/strike"+str(strike_value)+"_dip"+str(dip_value)+"/"
    ensure_dir(saved_path)

    ## check size of data
    datalen = len(farms_sliprate_data)

    ## slip
    plt.figure()
    plt.plot(benchmark_data[:,0],benchmark_data[:,1],'b-',label=benchmark_label)
    plt.plot(time[:datalen],farms_slip_data,'r--',label=farms_label)
    plt.title(benchmark_code+" slip time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ", fontsize=12)
    plt.legend()
    plt.xlabel("time (s)", fontsize=12)
    plt.ylabel("slip (m)", fontsize=12)
    plt.savefig(saved_path+"/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(benchmark_data[:,0],benchmark_data[:,2],'b-',label=benchmark_label)
    plt.plot(time[:datalen],farms_sliprate_data,'r--',label=farms_label)
    plt.title(benchmark_code+" slip rate time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ", fontsize=12)
    plt.legend()
    plt.xlabel("time (s)", fontsize=12)
    plt.ylabel("slip rate (m/s)", fontsize=12)
    plt.savefig(saved_path+"/sliprate.png")
    plt.show()   

    ## traction
    plt.figure()
    plt.plot(benchmark_data[:,0],benchmark_data[:,3],'b-',label=benchmark_label)
    plt.plot(time[:datalen],farms_traction_data,'r--',label=farms_label)
    plt.title(benchmark_code+"shear traction time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ",fontsize=12)
    plt.legend()
    plt.xlabel("time (s)", fontsize=12)
    plt.ylabel("shear traction (MPa)", fontsize=12)
    plt.savefig(saved_path+"/traction_x.png")
    plt.show()  

    ## normal traction
    plt.figure()
    plt.plot(benchmark_data[:,0],benchmark_data[:,-1],'b-',label=benchmark_label)
    plt.plot(time[:datalen],farms_traction_y,'r--',label=farms_label)
    plt.title(benchmark_code+"normal traction time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ",fontsize=12)
    plt.legend()
    plt.xlabel("time (s)", fontsize=12)
    plt.ylabel("normal traction (MPa)", fontsize=12)
    plt.savefig(saved_path+"/traction_y.png")
    plt.show()  

#strike,dip
#multiple nodes
given_coord_list = [[-2000 , 0, -7500],
                    [ 2000 , 0, -7500],
                    [ 5000 , 0, -7500],
                    [ 9000 , 0, -7500]]

# given_coord_list = [[   0  , -7500, 0]]

#num of files
num_of_file = np.shape(given_coord_list)[0]

#loop over files
for i in range(num_of_file):

    #get coords
    xcoord_i = given_coord_list[i][0] / 1000
    ycoord_i = given_coord_list[i][1] / 1000
    zcoord_i = given_coord_list[i][2] / 1000

    #file path
    benchmark_path = "./benchmark_data/benchmark_mainfault_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_slip_path = "./farms_data_elemental_qp/local_shear_jump_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_sliprate_path = "./farms_data_elemental_qp/local_shear_jump_rate_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_traction_path = "./farms_data_elemental_qp/local_shear_traction_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_tractiony_path = "./farms_data_elemental_qp/local_normal_traction_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"

    benchmark = np.loadtxt(benchmark_path, comments='#')
    farms_slip = np.loadtxt(farms_slip_path)
    farms_sliprate = np.loadtxt(farms_sliprate_path)
    farms_traction = np.loadtxt(farms_traction_path)
    farms_traction_y = np.loadtxt(farms_tractiony_path)

    print(farms_sliprate_path)

    ##plot
    plotfigure(benchmark,time,farms_sliprate,farms_slip,farms_traction,farms_traction_y,benchmark_label,farms_label,xcoord_i,zcoord_i,benchmark_code)