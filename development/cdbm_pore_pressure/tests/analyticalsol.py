#Analytical Solution of Fluid Injection
#Import packages
import numpy as np
from scipy.special import exp1
import matplotlib.pyplot as plt
import json

#---------------------------------#
#Define parameters
flux_q = 10; #kg/s
density_rho_0 = 1e3; #kg/m^3
permeability_k = 3e-12; #m^2
viscosity_eta = 0.4e-3; #Pa s
biotcoeff_alpha = 0.31; #-
undrained_nu_u = 0.3;  #-
shear_modulus_mu = 20e9; #Pa
drained_nu = 0.25; #-

#injection location
injection_x = 0.0
injection_y = 0.0

#geometry
nx = 100
ny = 100
xmin = -2500
xmax = 2500
ymin = -2500
ymax = 2500

#time
dt = 0.1
tmin = 0.0
tmax = 10.0
nt = int( ( tmax - tmin ) / dt )
#---------------------------------#

#generate data points
xdata = np.linspace(xmin,xmax,nx+1).round(0)
ydata = np.linspace(ymin,ymax,ny+1).round(0)
tdata = np.linspace(dt,tmax,nt).round(2)

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

#write
with open('analyticalsol.txt','w') as f:

    f.write('AXIS X')
    f.write('\n')
    for item in xdata_flat:
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
        # write each item on a new line
        f.write("%s " % item)
    f.write('\n')
    f.write('DATA')
    f.write('\n')
    for index_i in range(shapet):
        pdata_flat = arr_pressure_3d[:,:,index_i].flatten().tolist()
        for item in pdata_flat:
        # write each item on a new line
            f.write("%s\n" % item)


