import numpy as np
import matplotlib.pyplot as plt
import os

#benchmark code
benchmark_code = "TPV205"

#read benchmark data
benchmark_label = "benchmark-pylith-200m"

#read farms data
farms_label = "farms-200m"

#time farms 0.1s interval
time = np.linspace(0,12.0,121)

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
def plotfigure(benchmark_data, time, farms_sliprate_data, farms_slip_data, farms_traction_data, benchmark_label, farms_label, strike_value, dip_value, benchmark_code):

    ## check if dir exists, create one if not
    saved_path = "./outputs/strike"+str(strike_value)+"_dip"+str(dip_value)+"/"
    ensure_dir(saved_path)

    ## check size of data
    datalen = len(farms_sliprate_data)

    ## slip
    plt.figure()
    plt.plot(time[:datalen],-1.0*farms_slip_data,'g-',label=farms_label)
    plt.plot(benchmark_data[:,0],benchmark_data[:,1],'r-',label=benchmark_label)
    plt.title(benchmark_code+" slip time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ")
    plt.legend()
    plt.xlabel("time (s)")
    plt.ylabel("slip (m)")
    plt.savefig(saved_path+"/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time[:datalen],-1.0*farms_sliprate_data,'g-',label=farms_label)
    plt.plot(benchmark_data[:,0],benchmark_data[:,2],'r-',label=benchmark_label)
    plt.title(benchmark_code+" slip rate time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ")
    plt.legend()
    plt.xlabel("time (s)")
    plt.ylabel("slip rate (m/s)")
    plt.savefig(saved_path+"/sliprate.png")
    plt.show()   

    ## traction
    plt.figure()
    plt.plot(time[:datalen],-1.0*farms_traction_data/1e6+benchmark_data[1,3],'g-',label=farms_label)
    plt.plot(benchmark_data[:,0],benchmark_data[:,3],'r-',label=benchmark_label)
    plt.title(benchmark_code+" traction time history at strike "+str(strike_value)+"km and at dip "+str(dip_value)+"km ")
    plt.legend()
    plt.xlabel("time (s)")
    plt.ylabel("traction (MPa)")
    plt.savefig(saved_path+"/traction.png")
    plt.show()  

#strike,dip
#multiple nodes
given_coord_list = [[ 4500 , 0, -7500],
                    [-7500 , 0, -7500],
                    [-12000, 0, -7500],
                    [   0  , 0, -7500],              
                    [4500  , 0, -7500],
                    [7500  , 0, -7500],
                    [12000 , 0, -7500]]

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
    benchmark_path = "./benchmark_data/pylith_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_slip_path = "./farms_data_elemental/jump_x_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_sliprate_path = "./farms_data_elemental/jump_x_rate_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"
    farms_traction_path = "./farms_data_elemental/traction_x_strike"+str(xcoord_i)+"_dip"+str(zcoord_i)+".txt"

    benchmark = np.loadtxt(benchmark_path)
    farms_slip = np.loadtxt(farms_slip_path)
    farms_sliprate = np.loadtxt(farms_sliprate_path)
    farms_traction = np.loadtxt(farms_traction_path)

    print(farms_sliprate_path)

    ##plot
    plotfigure(benchmark,time,farms_sliprate,farms_slip,farms_traction,benchmark_label,farms_label,xcoord_i,zcoord_i,benchmark_code)