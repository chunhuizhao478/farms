import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-1/x)

# Focus on x > 0, in a range where the function's gradient is visible
x = np.linspace(0.001, 50, 500)  # adjust as desired
y = f(x)

plt.figure(figsize=(7, 5))
plt.plot(x, y, label=r'$y = e^{-1/x}, \; x > 0$')
plt.xlabel('T', fontsize=14)
plt.ylabel(r'$y = e^{-1/T}$', fontsize=14)
plt.title('Plot of $y = e^{-1/T}$ for $T > 0$', fontsize=16)
plt.grid(True)
plt.legend()
plt.show()