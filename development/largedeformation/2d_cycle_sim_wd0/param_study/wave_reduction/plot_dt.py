import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("dynamic_solve_damping_out.csv",skiprows=1,delimiter=',')

plt.figure()
plt.plot(data[90:,0],data[90:,1],'r-*')
plt.title("Time History of Time Step Evolution")
plt.ylabel("Time Step (s)")
plt.xlabel("Time (s)")
plt.show()