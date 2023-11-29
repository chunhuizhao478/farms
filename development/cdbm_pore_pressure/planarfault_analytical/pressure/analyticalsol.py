#Analytical Solution of Fluid Injection
#Reference: RUDNICKI, FLUID MASS SOURCES AND POINT FORCES IN LINEAR ELASTIC DIFFUSIVE SOLIDS 
#Import packages
import numpy as np
from scipy.special import exp1
import matplotlib.pyplot as plt
import json

#---------------------------------#
#Define parameters
flux_q = 1e4; #kg/s
density_rho_0 = 1e3; #kg/m^3
permeability_k = 3e-12; #m^2
viscosity_eta = 0.4e-3; #Pa s
biotcoeff_alpha = 0.31; #-
undrained_nu_u = 0.3;  #-
shear_modulus_mu = 32.04e9; #Pa
drained_nu = 0.25; #-

#injection location
injection_x = 0.0
injection_y = 0.0

#geometry
nx = 480
ny = 160
xmin = -6000
xmax = 6000
ymin = -2000
ymax = 2000

#time
dt = 5e-4
tmin = 0.0
tmax = dt * 10000
nt = int( ( tmax - tmin ) / dt )
#---------------------------------#

#generate data points
xdata = np.linspace(xmin,xmax,nx+1).round(0)
ydata = np.linspace(ymin,ymax,ny+1).round(0)
tdata = np.linspace(dt,tmax,nt).round(5) #!

# print(tdata)
# exit()

xv, yv = np.meshgrid(xdata,ydata)

shapexv0 = np.shape(xv)[0]
shapexv1 = np.shape(xv)[1]
shapet = np.shape(tdata)[0]

# print(shapet)
# exit()

#compute R (w.r.t the injection location)
arr_injection_x = np.ones(np.shape(xv)) * injection_x
arr_injection_y = np.ones(np.shape(yv)) * injection_y

arr_diff_x = xv - arr_injection_x
arr_diff_y = yv - arr_injection_y

arr_R = np.sqrt(arr_diff_x*arr_diff_x+arr_diff_y*arr_diff_y)

# print(np.max(arr_R),np.sqrt(2)*2500)
# print(np.min(arr_R))

#undrained lame constant
drained_lambda   = 2 * shear_modulus_mu *     drained_nu / ( 1 - 2 *     drained_nu )
undrained_lambda = 2 * shear_modulus_mu * undrained_nu_u / ( 1 - 2 * undrained_nu_u )

#hydraulic diffusivity
c = ( permeability_k * ( undrained_lambda - drained_lambda ) * ( drained_lambda + 2 * shear_modulus_mu ) ) / ( viscosity_eta * biotcoeff_alpha * biotcoeff_alpha * ( undrained_lambda + 2 * shear_modulus_mu ) )

#compute pressure
arr_pressure_3d = np.zeros((shapexv0,shapexv1,shapet))

#plot figures
# plt.figure()
index = 0
for t_i in tdata[:]:

    print("Time: ", t_i)

    arr_pressure_3d[:,:,index] = ( flux_q * viscosity_eta ) / ( 4 * np.pi * density_rho_0 * permeability_k ) * exp1( arr_R * arr_R / ( 4 * c * t_i ) )

    arr_pressure_3d[:,:,index] = arr_pressure_3d[:,:,index].round(4)

    index += 1
    # plt.imshow(arr_pressure_3d[:,:,t_i], cmap='cool', interpolation='nearest')
    # plt.colorbar()
    # plt.show()
    # plt.close()

#save file into txt

#rearrange xdata ydata
xdata_flat = xdata.flatten().tolist()
ydata_flat = ydata.flatten().tolist()
tdata_flat = tdata.tolist()

#set a threshold to inf in the middle point (0.9*sigma_yy)

#write
with open('analyticalsol.txt','w') as f:

    f.write('AXIS X')
    f.write('\n')
    for item in xdata_flat:

        print("save AXIS X ...")

        # write each item on a new line
        f.write("%s " % item)
    f.write('\n')
    f.write('AXIS Y')
    f.write('\n')
    for item in ydata_flat:
        # write each item on a new line
        f.write("%s " % item)
    f.write('\n')
    f.write('AXIS T')
    f.write('\n')
    for item in tdata_flat:

        print("save AXIS Y ...")

        # write each item on a new line
        f.write("%s " % item)
    f.write('\n')
    f.write('DATA')
    f.write('\n')
    count = 0
    for index_i in range(shapet):
        pdata_flat = arr_pressure_3d[:,:,index_i].flatten().tolist()

        count += 1
        print("save DATA ... ", round(count/shapet*100, 3), " %")
        for item in pdata_flat:

            #check if inf, set a threshold
            boolreturn = np.isinf(item)

            if boolreturn:

                item = 0.9 * 120e6

        # write each item on a new line
            f.write("%s\n" % item)



