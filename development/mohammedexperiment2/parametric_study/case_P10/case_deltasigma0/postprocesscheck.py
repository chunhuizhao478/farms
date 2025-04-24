import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pd_data = pd.read_csv("test_me_maintb_csv_main_fault_1440.csv")

# material properties
mu_s = 0.7  # static friction coefficient
mu_d = 0.1  # dynamic friction coefficient

# Extracting the relevant columns
T1 = pd_data['T1_aux']
T2 = pd_data['T2_aux']
jump_x_aux = pd_data['jump_x_aux']
jump_y_aux = pd_data['jump_y_aux']
jump_x_rate_aux = pd_data['jump_x_rate_aux']
jump_y_rate_aux = pd_data['jump_y_rate_aux']
x = pd_data['x']
y = pd_data['y']

# get coordinate of the main fault
arclen_x = x - x[0]
arclen_y = y - y[0]
arclen = np.sqrt(arclen_x**2 + arclen_y**2)

# get local jump x for main fault
# local_jump = jump_x_aux * cos(29 degree) + jump_y_aux * sin(29 degree)
local_jump_x = jump_x_aux * np.cos(np.radians(29)) + jump_y_aux * np.sin(np.radians(29))
local_jump_rate_x = jump_x_rate_aux * np.cos(np.radians(29)) + jump_y_rate_aux * np.sin(np.radians(29))

#plot shear traction
plt.figure(figsize=(10, 5))
plt.plot(arclen, -T1, 'r-*' , label='shear traction (MPa)')
plt.plot(arclen, -T2 * mu_s, 'b-*', label='peak shear strength (MPa)')
plt.plot(arclen, -T2 * mu_d, 'k-*', label='residual shear strength (MPa)')
plt.legend()
plt.xlabel('arc length (m)')
plt.ylabel('shear traction (MPa)')

plt.figure(figsize=(10, 5))
plt.plot(arclen, local_jump_x * 1e3, 'r-*', label='local tangent jump')
plt.legend()
plt.xlabel('arc length (mm)')
plt.ylabel('jump (m)')
plt.title('tangent jump')

plt.figure(figsize=(10, 5))
plt.plot(arclen, local_jump_rate_x, 'r-*', label='local tangent jump rate')
plt.legend()
plt.xlabel('arc length (m)')
plt.ylabel('tangent jump rate (m/s)')
plt.title('tangent jump rate')

plt.show()