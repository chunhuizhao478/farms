import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("dynamic_solve_csv.csv",skiprows=1,delimiter=',')

plt.figure()
plt.loglog(data1[82:,0],data1[82:,1],'r-*',label="Quadratic, w damper, QuasiDyna")
plt.title("Time History of Time Step Evolution")
plt.ylabel("Time Step (s)")
plt.xlabel("Time (s)")
plt.legend()
plt.show()