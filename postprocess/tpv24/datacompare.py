import numpy as np
import matplotlib.pyplot as plt

#read benchmark data
strikem8dip10 = np.loadtxt("./benchmark_data/eqdynastrikem8dip10.txt")

#
strikem8dip10_2 = np.loadtxt("./benchmark_data/seisolstrikem8dip10.txt")

strikem8dip5 = np.loadtxt("./benchmark_data/eqdynastrikem8dip5.txt")

strikem8dip0 = np.loadtxt("./benchmark_data/eqdynastrikem8dip0.txt")

strikem2dip10 = np.loadtxt("./benchmark_data/eqdynastrikem2dip10.txt")

#read farms data
sliprate_strikem8dip10 = np.loadtxt("./farms_data/sliprate_strikem8_dip10.txt",skiprows=1)
slip_strikem8dip10 = np.loadtxt("./farms_data/slip_strikem8_dip10.txt",skiprows=1)

#
sliprate_strikem8dip10_new = np.loadtxt("./farms_data/jun11/sliprate.txt",skiprows=1)
slip_strikem8dip10_new = np.loadtxt("./farms_data/jun11/slip.txt",skiprows=1)

#
sliprate_strikem8dip10_new2 = np.loadtxt("./farms_data/jun13/sliprate.txt",skiprows=1)
slip_strikem8dip10_new2 = np.loadtxt("./farms_data/jun13/slip.txt",skiprows=1)

sliprate_strikem8dip5 = np.loadtxt("./farms_data/sliprate_strikem8_dip5.txt",skiprows=1)
slip_strikem8dip5 = np.loadtxt("./farms_data/slip_strikem8_dip5.txt",skiprows=1)

sliprate_strikem8dip0 = np.loadtxt("./farms_data/sliprate_strikem8_dip0.txt",skiprows=1)
slip_strikem8dip0 = np.loadtxt("./farms_data/slip_strikem8_dip0.txt",skiprows=1)

sliprate_strikem2dip10 = np.loadtxt("./farms_data/sliprate_strikem2_dip10.txt",skiprows=1)
slip_strikem2dip10 = np.loadtxt("./farms_data/slip_strikem2_dip10.txt",skiprows=1)

#time farms
time = np.linspace(0,12.0,121)

plot_strikem8dip10 = True
plot_strikem8dip5 = False
plot_strikem8dip0 = False
plot_strikem2dip10 = False

# plot
# strikem8dip10
if plot_strikem8dip10:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem8dip10,'k-',label="farms-old")
    plt.plot(time[:88],slip_strikem8dip10_new,'b-',label="farms-high damp")
    plt.plot(time[:76],slip_strikem8dip10_new2,'m-',label="farms-high damp-high mus")
    plt.plot(strikem8dip10[:,0],strikem8dip10[:,1],'r-',label="benchmark-eqdyna")
    plt.plot(strikem8dip10_2[:,0],2*strikem8dip10_2[:,1],'r-',label="benchmark-seisol")
    plt.title("TPV24 slip time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip10/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem8dip10,'k-',label="farms-old")
    plt.plot(time[:88],sliprate_strikem8dip10_new,'b-',label="farms-high damp")
    plt.plot(time[:76],sliprate_strikem8dip10_new2,'m-',label="farms-high damp-high mus")
    plt.plot(strikem8dip10[:,0],strikem8dip10[:,2],'r-',label="benchmark-eqdyna")
    plt.plot(strikem8dip10_2[:,0],strikem8dip10_2[:,2],'r-',label="benchmark-seisol")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip10/sliprate.png")
    plt.show()

# strikem8dip5
if plot_strikem8dip5:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem8dip5,'k-',label="farms")
    plt.plot(strikem8dip5[:,0],strikem8dip5[:,1],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip time history at strike -8km and at dip 5km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip5/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem8dip5,'k-',label="farms")
    plt.plot(strikem8dip5[:,0],strikem8dip5[:,2],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 5km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip5/sliprate.png")
    plt.show()

# strikem8dip5
if plot_strikem8dip0:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem8dip0,'k-',label="farms")
    plt.plot(strikem8dip0[:,0],strikem8dip0[:,1],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip time history at strike -8km and at dip 0km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip0/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem8dip0,'k-',label="farms")
    plt.plot(strikem8dip0[:,0],strikem8dip0[:,2],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 0km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip0/sliprate.png")
    plt.show()

# strikem2dip10
if plot_strikem2dip10:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem2dip10,'k-',label="farms")
    plt.plot(strikem2dip10[:,0],strikem2dip10[:,1],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip time history at strike -2km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem2dip10/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem2dip10,'k-',label="farms")
    plt.plot(strikem2dip10[:,0],strikem2dip10[:,2],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip rate time history at strike -2km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem2dip10/sliprate.png")
    plt.show()