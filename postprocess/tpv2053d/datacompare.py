import numpy as np
import matplotlib.pyplot as plt

#read benchmark data
strikem0dip7dot5 = np.loadtxt("./benchmark_data/eqdyna_strike0_dip7dot5.txt")

#read farms data
sliprate_strike0dip7dot5 = np.loadtxt("./farms_data/sliprate_strike0_dip7dot5.txt",skiprows=1)
slip_strike0dip7dot5 = np.loadtxt("./farms_data/slip_strike0_dip7dot5.txt",skiprows=1)

#time farms
time = np.linspace(0,12.0,121)

plot_strike0_dip7dot5 = True
# plot_strikem8dip5 = True
# plot_strikem8dip0 = True
# plot_strikem2dip10 = True

# plot
# strike0_dip7dot5
if plot_strike0_dip7dot5:

    ## slip
    plt.figure()
    plt.plot(time[:32],slip_strike0dip7dot5,'g-',label="farms-high damp-high mus-2")
    plt.plot(strikem0dip7dot5[:,0],strikem0dip7dot5[:,1],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strike0_dip7dot5/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time[:32],sliprate_strike0dip7dot5,'g-',label="farms-high damp-high mus-2")
    plt.plot(strikem0dip7dot5[:,0],strikem0dip7dot5[:,2],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strike0_dip7dot5/sliprate.png")
    plt.show()

    ## shear stress
    # plt.figure()
    # plt.plot(time[:73],shearsts_strikem8dip10_new4,'k-',label="farms-high damp-high mus-2")
    # plt.plot(strikem8dip10[:,0],strikem8dip10[:,3],'r-',label="benchmark-eqdyna")
    # plt.title("TPV24 shear stress time history at strike -8km and at dip 10km ")
    # plt.legend()
    # plt.ylabel("shear stress (MPa)")
    # plt.xlabel("time (s)")
    # plt.savefig("./outputs/strikem8dip10/shearstress.png")
    # plt.show()