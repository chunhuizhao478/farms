import meshio_run_tria as func

#Generate properties files 
#Compression 
P = [-10e6, -15e6, -20e6]; 
P_label = ['10','15','20'] #MPa
nu = 0.37 #poisson ratio, the lateral confining pressure: nu * P
delta_normalsts = [-0e6, -2e6, -4e6, -6e6]; 
delta_normalsts_label = ['0', '2', '4', '6'] #MPa
meshfilepath = './mesh/me.msh'

# P = [-10e6]; 
# P_label = ['10'] #MPa
# nu = 0.37 #poisson ratio, the lateral confining pressure: nu * P
# delta_normalsts = [-0e6]; 
# delta_normalsts_label = ['0'] #MPa
# meshfilepath = './mesh/me.msh'

#---------------------------------------------------------------#
#Point at distance 0.015 from P9: (0.134155, 0.129788, 0.000000)#

branch_end_point = [0.134155, 0.129788, 0.000000] #x,y,z
#---------------------------------------------------------------#

for i in range(len(P)):
    for j in range(len(delta_normalsts)):
       print(f'Running for P = {P[i]} and delta_normalsts = {delta_normalsts[j]}')
       func.main_func(sigmayy= P[i], 
                      nu=nu, 
                      added_normalsts=delta_normalsts[j],
                      sigmayy_label=P_label[i],
                      added_normalsts_label=delta_normalsts_label[j],
                      mesh_path=meshfilepath,
                      branch_end_point=branch_end_point)

