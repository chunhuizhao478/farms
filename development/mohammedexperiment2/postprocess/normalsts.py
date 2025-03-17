import numpy as np
import matplotlib.pyplot as plt

normalsts = [3, 6, 9, 12]
time = [0.000189342, 0.000198387, 0.000209844, 0.000219492]

plt.figure(figsize=(7, 5))
plt.plot(normalsts, time, label='Time of Rupture Nucleation vs Normal stresses')
plt.xlabel('Normal stresses', fontsize=14)
plt.ylabel('Time', fontsize=14)
plt.title('Time vs Normal stresses', fontsize=16)
plt.grid(True)
plt.legend()
plt.show()
