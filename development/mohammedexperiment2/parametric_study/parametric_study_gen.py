import meshio_run_tria as func

#Generate properties files 
#Compression 
P = [-10e6, -20e6]; 
P_label = ['10', '20'] #MPa
nu = 0.37 #poisson ratio, the lateral confining pressure: nu * P
delta_normalsts = [-0e6, -2e6, -4e6, -6e6, -8e6, -10e6]; 
delta_normalsts_label = ['0', '2', '4', '6', '8', '10'] #MPa
meshfilepath = './mesh/me.msh'

for i in range(len(P)):
    for j in range(len(delta_normalsts)):
       print(f'Running for P = {P[i]} and delta_normalsts = {delta_normalsts[j]}')
       func.main_func(sigmayy= P[i], 
                      nu=nu, 
                      added_normalsts=delta_normalsts[j],
                      sigmayy_label=P_label[i],
                      added_normalsts_label=delta_normalsts_label[j],
                      mesh_path=meshfilepath)

