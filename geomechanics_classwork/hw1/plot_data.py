import numpy as np
import matplotlib.pyplot as plt

#read data from file
barnett_density_data = np.loadtxt('barnett_density_data.txt', skiprows=1)

plt.figure(figsize=(10, 6))
plt.plot(barnett_density_data[:, 1], barnett_density_data[:, 0], '.-', label='Barnett Shale Density')
plt.xlabel("Density (g/cm³)")
plt.ylabel("Depth (ft)")
plt.gca().invert_yaxis()
plt.title("Barnett Shale Density vs Depth")
plt.grid()
plt.legend()
plt.savefig('barnett_density_plot.png', dpi=300)
plt.show()

#read data from file
GOM_offshore_data = np.loadtxt('GOM_offshore_data.txt', skiprows=1)
plt.figure(figsize=(10, 6))
plt.plot(GOM_offshore_data[:, 1], GOM_offshore_data[:, 0], '.-', label='GOM Offshore Density')
plt.xlabel("Density (g/cm³)")
plt.ylabel("Depth (ft)")
plt.gca().invert_yaxis()
plt.title("GOM Offshore Density vs Depth")
plt.grid()
plt.legend()
plt.savefig('GOM_offshore_density_plot.png', dpi=300)
plt.show()