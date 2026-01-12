import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pd_data = pd.read_csv("tpv142D_tria_1_csv_main_fault_0760.csv")

local_shear_traction = pd_data['local_shear_traction']
fractal_stress = pd_data['fractal_shear_stress_aux']

plt.figure()
plt.plot(local_shear_traction, label='local_shear_traction')
# plt.plot(fractal_stress, label='local_shear_traction')
plt.show() 