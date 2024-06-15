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

#jun11
sliprate_strikem8dip10_new = np.loadtxt("./farms_data/jun11/sliprate.txt",skiprows=1)
slip_strikem8dip10_new = np.loadtxt("./farms_data/jun11/slip.txt",skiprows=1)

#jun13
sliprate_strikem8dip10_new2 = np.loadtxt("./farms_data/jun13/sliprate.txt",skiprows=1)
slip_strikem8dip10_new2 = np.loadtxt("./farms_data/jun13/slip.txt",skiprows=1)

#jun14
sliprate_strikem8dip10_new3 = np.loadtxt("./farms_data/jun14/slipratestrikem8dip10.txt",skiprows=1)
slip_strikem8dip10_new3 = np.loadtxt("./farms_data/jun14/slipstrikem8dip10.txt",skiprows=1)
sliprate_strikem8dip5_new3 = np.loadtxt("./farms_data/jun14/slipratestrikem8dip5.txt",skiprows=1)
slip_strikem8dip5_new3 = np.loadtxt("./farms_data/jun14/slipstrikem8dip5.txt",skiprows=1)
sliprate_strikem8dip0_new3 = np.loadtxt("./farms_data/jun14/slipratestrikem8dip0.txt",skiprows=1)
slip_strikem8dip0_new3 = np.loadtxt("./farms_data/jun14/slipstrikem8dip0.txt",skiprows=1)
sliprate_strikem2dip10_new3 = np.loadtxt("./farms_data/jun14/slipratestrikem2dip10.txt",skiprows=1)
slip_strikem2dip10_new3 = np.loadtxt("./farms_data/jun14/slipstrikem2dip10.txt",skiprows=1)

#jun15
sliprate_strikem8dip10_new4 = np.loadtxt("./farms_data/jun15/slipratestrikem8dip10.txt",skiprows=1)
slip_strikem8dip10_new4 = np.loadtxt("./farms_data/jun15/slipstrikem8dip10.txt",skiprows=1)
shearsts_strikem8dip10_new4 = np.loadtxt("./farms_data/jun15/shearstressm8dip10.txt",skiprows=1)
shearsts_strikem8dip5_new4 = np.loadtxt("./farms_data/jun15/shearstressm8dip5.txt",skiprows=1)
shearsts_strikem8dip0_new4 = np.loadtxt("./farms_data/jun15/shearstressm8dip0.txt",skiprows=1)
shearsts_strikem2dip10_new4 = np.loadtxt("./farms_data/jun15/shearstressm2dip10.txt",skiprows=1)

sliprate_strikem8dip5 = np.loadtxt("./farms_data/sliprate_strikem8_dip5.txt",skiprows=1)
slip_strikem8dip5 = np.loadtxt("./farms_data/slip_strikem8_dip5.txt",skiprows=1)

sliprate_strikem8dip0 = np.loadtxt("./farms_data/sliprate_strikem8_dip0.txt",skiprows=1)
slip_strikem8dip0 = np.loadtxt("./farms_data/slip_strikem8_dip0.txt",skiprows=1)

sliprate_strikem2dip10 = np.loadtxt("./farms_data/sliprate_strikem2_dip10.txt",skiprows=1)
slip_strikem2dip10 = np.loadtxt("./farms_data/slip_strikem2_dip10.txt",skiprows=1)

#time farms
time = np.linspace(0,12.0,121)

plot_strikem8dip10 = True
plot_strikem8dip5 = True
plot_strikem8dip0 = True
plot_strikem2dip10 = True

# plot
# strikem8dip10
if plot_strikem8dip10:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem8dip10,'k-',label="farms-old")
    plt.plot(time[:76],slip_strikem8dip10_new3,'g-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip10[:,0],strikem8dip10[:,1],'r-',label="benchmark-eqdyna")
    plt.plot(strikem8dip10_2[:,0],2*strikem8dip10_2[:,1],'r-',label="benchmark-seisol")
    plt.title("TPV24 slip time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip10/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem8dip10,'k-',label="farms-old")
    plt.plot(time[:76],sliprate_strikem8dip10_new3,'g-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip10[:,0],strikem8dip10[:,2],'r-',label="benchmark-eqdyna")
    plt.plot(strikem8dip10_2[:,0],strikem8dip10_2[:,2],'r-',label="benchmark-seisol")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip10/sliprate.png")
    plt.show()

    ## shear stress
    plt.figure()
    plt.plot(time[:73],shearsts_strikem8dip10_new4,'k-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip10[:,0],strikem8dip10[:,3],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 shear stress time history at strike -8km and at dip 10km ")
    plt.legend()
    plt.ylabel("shear stress (MPa)")
    plt.xlabel("time (s)")
    plt.savefig("./outputs/strikem8dip10/shearstress.png")
    plt.show()

# strikem8dip5
if plot_strikem8dip5:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem8dip5,'k-',label="farms")
    plt.plot(time[:76],slip_strikem8dip5_new3,'g-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip5[:,0],strikem8dip5[:,1],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip time history at strike -8km and at dip 5km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip5/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem8dip5,'k-',label="farms")
    plt.plot(time[:76],sliprate_strikem8dip5_new3,'g-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip5[:,0],strikem8dip5[:,2],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 5km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip5/sliprate.png")
    plt.show()

    ## shear stress
    plt.figure()
    plt.plot(time[:73],shearsts_strikem8dip5_new4,'k-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip5[:,0],strikem8dip5[:,3],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 shear stress time history at strike -8km and at dip 5km ")
    plt.legend()
    plt.ylabel("shear stress (MPa)")
    plt.xlabel("time (s)")
    plt.savefig("./outputs/strikem8dip5/shearstress.png")
    plt.show()

# strikem8dip0
if plot_strikem8dip0:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem8dip0,'k-',label="farms")
    plt.plot(strikem8dip0[:,0],strikem8dip0[:,1],'r-',label="benchmark-eqdyna")
    plt.plot(time[:76],slip_strikem8dip0_new3,'g-',label="farms-high damp-high mus-2")
    plt.title("TPV24 slip time history at strike -8km and at dip 0km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip0/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem8dip0,'k-',label="farms")
    plt.plot(strikem8dip0[:,0],strikem8dip0[:,2],'r-',label="benchmark-eqdyna")
    plt.plot(time[:76],sliprate_strikem8dip0_new3,'g-',label="farms-high damp-high mus-2")
    plt.title("TPV24 slip rate time history at strike -8km and at dip 0km ")
    plt.legend()
    plt.savefig("./outputs/strikem8dip0/sliprate.png")
    plt.show()

    ## shear stress
    plt.figure()
    plt.plot(time[:73],shearsts_strikem8dip0_new4,'k-',label="farms-high damp-high mus-2")
    plt.plot(strikem8dip0[:,0],strikem8dip0[:,3],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 shear stress time history at strike -8km and at dip 0km ")
    plt.legend()
    plt.ylabel("shear stress (MPa)")
    plt.xlabel("time (s)")
    plt.savefig("./outputs/strikem8dip0/shearstress.png")
    plt.show()

# strikem2dip10
if plot_strikem2dip10:

    ## slip
    plt.figure()
    plt.plot(time,slip_strikem2dip10,'k-',label="farms")
    plt.plot(strikem2dip10[:,0],strikem2dip10[:,1],'r-',label="benchmark-eqdyna")
    plt.plot(time[:76],slip_strikem2dip10_new3,'g-',label="farms-high damp-high mus-2")
    plt.title("TPV24 slip time history at strike -2km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem2dip10/slip.png")
    plt.show()

    ## slip rate
    plt.figure()
    plt.plot(time,sliprate_strikem2dip10,'k-',label="farms")
    plt.plot(strikem2dip10[:,0],strikem2dip10[:,2],'r-',label="benchmark-eqdyna")
    plt.plot(time[:76],sliprate_strikem2dip10_new3,'g-',label="farms-high damp-high mus-2")
    plt.title("TPV24 slip rate time history at strike -2km and at dip 10km ")
    plt.legend()
    plt.savefig("./outputs/strikem2dip10/sliprate.png")
    plt.show()

    ## shear stress
    plt.figure()
    plt.plot(time[:73],shearsts_strikem2dip10_new4,'k-',label="farms-high damp-high mus-2")
    plt.plot(strikem2dip10[:,0],strikem2dip10[:,3],'r-',label="benchmark-eqdyna")
    plt.title("TPV24 shear stress time history at strike -2km and at dip 0km ")
    plt.legend()
    plt.ylabel("shear stress (MPa)")
    plt.xlabel("time (s)")
    plt.savefig("./outputs/strikem2dip10/shearstress.png")
    plt.show()