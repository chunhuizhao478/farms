import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./segment_info_array.txt")

ptrA_xcoord = data[:,0]
ptrA_ycoord = data[:,1]
ptrB_xcoord = data[:,2]
ptrB_ycoord = data[:,3]

local_shear_sts = data[:,4]
local_normal_sts = data[:,5]

mus = 0.677
mud = 0.4

pressure_needed = (abs(local_shear_sts) / mus - abs(local_normal_sts))

mu_values = data[:,-2]
S_values = data[:,-1]

plt.figure()
for i in range(len(ptrA_xcoord)):
    if S_values[i] < 0:
        if  mu_values[i] > 0.677:
            plt.plot([ptrA_xcoord[i], ptrB_xcoord[i]],[ptrA_ycoord[i], ptrB_ycoord[i]],'m-')
        else:
            plt.plot([ptrA_xcoord[i], ptrB_xcoord[i]],[ptrA_ycoord[i], ptrB_ycoord[i]],'k-')
    else:
        plt.plot([ptrA_xcoord[i], ptrB_xcoord[i]],[ptrA_ycoord[i], ptrB_ycoord[i]],'r-')
    # plt.text(0.5*(ptrA_xcoord[i]+ptrB_xcoord[i]), 0.5*(ptrA_ycoord[i]+ptrB_ycoord[i]), f'{S_values[i]}', ha='right')
    plt.text(0.5*(ptrA_xcoord[i]+ptrB_xcoord[i]), 0.5*(ptrA_ycoord[i]+ptrB_ycoord[i]), f'{pressure_needed[i]:.2f}', ha='right')

plt.title("S value distribution")
plt.xlabel("coordinate x (m)")
plt.ylabel("coordinate y (m)")
plt.show()