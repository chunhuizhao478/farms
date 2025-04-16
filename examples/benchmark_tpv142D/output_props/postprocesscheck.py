import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pd_data = pd.read_csv("tpv142D_tria_csv_main_fault_0500.csv")

local_shear_traction = pd_data['local_shear_traction']

plt.figure()
plt.plot(local_shear_traction, label='local_shear_traction')
plt.show()