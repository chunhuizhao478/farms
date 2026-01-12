import pandas as pd

pd_data1 = pd.read_csv("main_fault_fractal_stress_0.csv")
pd_data2 = pd.read_csv("main_fault_fractal_stress_1.csv")
pd_data3 = pd.read_csv("main_fault_fractal_stress_2.csv")
        
# Optionally, visualize the profile (requires matplotlib)
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.plot(pd_data1['shear_stress'])
plt.plot(pd_data2['shear_stress'])
plt.plot(pd_data3['shear_stress'])
plt.title("Fractal Shear Stress Profile")
plt.xlabel("Index")
plt.ylabel("Shear Stress")
plt.show()