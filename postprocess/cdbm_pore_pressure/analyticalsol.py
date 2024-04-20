#Analytical Solution of Fluid Injection
#Reference: RUDNICKI, FLUID MASS SOURCES AND POINT FORCES IN LINEAR ELASTIC DIFFUSIVE SOLIDS 
#Import packages
import numpy as np
from scipy.special import exp1, erfc
import matplotlib.pyplot as plt
import json

#---------------------------------#
#Define parameters
flux_q = 10; #kg/s
density_rho_0 = 1e3; #kg/m^3
permeability_k = 3e-16; #m^2
viscosity_eta = 0.4e-3; #Pa s
biotcoeff_alpha = 0.31; #-
undrained_nu_u = 0.3;  #-
shear_modulus_mu = 32.04e9; #Pa
drained_nu = 0.25; #-

#injection location
injection_x = 0.0
injection_y = 0.0

#geometry
nx = 2
ny = 2
xmin = -1e-1
xmax = 1e-1
ymin = -1e-1
ymax = 1e-1

#time
dt = 0.5 * 24 * 60 * 60
tmin = 0.0
tmax = 1 * 24 * 60 * 60
nt = int( ( tmax - tmin ) / dt )
#---------------------------------#

tdata = np.linspace(dt,tmax,nt).round(2)

shapexv0 = 1
shapexv1 = 1
shapet = np.shape(tdata)[0]

arr_diff_x = 10
arr_diff_y = 10

arr_R = np.sqrt(arr_diff_x*arr_diff_x+arr_diff_y*arr_diff_y)

# print(np.max(arr_R),np.sqrt(2)*2500)
# print(np.min(arr_R))

#undrained lame constant
drained_lambda   = 2 * shear_modulus_mu *     drained_nu / ( 1 - 2 *     drained_nu )
undrained_lambda = 2 * shear_modulus_mu * undrained_nu_u / ( 1 - 2 * undrained_nu_u )

#hydraulic diffusivity
c = ( permeability_k * ( undrained_lambda - drained_lambda ) * ( drained_lambda + 2 * shear_modulus_mu ) ) / ( viscosity_eta * biotcoeff_alpha * biotcoeff_alpha * ( undrained_lambda + 2 * shear_modulus_mu ) )

print(c)
# exit()

#compute pressure
arr_pressure_3d = np.zeros((shapet,1))

#plot figures
# plt.figure()
index = 0
for t_i in tdata[:]:

    # arr_pressure_3d[index,0] = ( flux_q * viscosity_eta ) / ( 4 * np.pi * density_rho_0 * permeability_k ) * exp1( arr_R * arr_R / ( 4 * c * t_i ) )

    # 2D
    # arr_pressure_3d[index,0] = arr_pressure_3d[index,0].round(4)

    #3D
    xi = arr_R / np.sqrt(c*t_i)
    arr_pressure_3d[index,0] = ( flux_q * viscosity_eta ) / ( 4 * np.pi * density_rho_0 * arr_R * permeability_k ) * erfc( 0.5 * xi )

    index += 1
    # plt.imshow(arr_pressure_3d[:,:,t_i], cmap='cool', interpolation='nearest')
    # plt.colorbar()
    # plt.show()
    # plt.close()

#save file into txt

#rearrange xdata ydata
tdata_flat = tdata.tolist()
arr_pressure_3d = arr_pressure_3d.tolist()

np.savetxt("./tdata_flat.txt",
                tdata_flat,
                newline=" ")
np.savetxt('./pressure.txt',
                arr_pressure_3d,
                newline=" ")

# #write
# with open('analyticalsol.txt','w') as f:

#     f.write('AXIS X')
#     f.write('\n')
#     for item in xdata_flat:
#         # write each item on a new line
#         f.write("%s " % item)
#     f.write('\n')
#     f.write('AXIS Y')
#     f.write('\n')
#     for item in ydata_flat:
#         # write each item on a new line
#         f.write("%s " % item)
#     f.write('\n')
#     f.write('AXIS T')
#     f.write('\n')
#     for item in tdata_flat:
#         # write each item on a new line
#         f.write("%s " % item)
#     f.write('\n')
#     f.write('DATA')
#     f.write('\n')
#     for index_i in range(shapet):
#         pdata_flat = arr_pressure_3d[:,:,index_i].flatten().tolist()
#         for item in pdata_flat:
#         # write each item on a new line
#             f.write("%s\n" % item)



