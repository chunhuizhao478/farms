import pandas as pd
import matplotlib.pyplot as plt

pd_moose = pd.read_csv("../../fractal_stress/tpv142D_tria_csv_main_fault_0020.csv")
pd_moose_data = pd_moose[['fractal_shear_stress_aux']]

pd_precomputed = pd.read_csv("main_fault_fractal_stress.csv")
pd_precomputed_data = pd_precomputed[['shear_stress']]

# Plotting the data
plt.figure(figsize=(10, 6))
plt.plot(pd_moose_data, label='Moose Data', color='blue')
plt.plot(pd_precomputed_data/1e6, label='Precomputed Data', color='orange')
plt.title('Comparison of Fractal Shear Stress')
plt.xlabel('Index')
plt.ylabel('Fractal Shear Stress')  
plt.legend()
plt.grid()
plt.show()

#verify the data :)