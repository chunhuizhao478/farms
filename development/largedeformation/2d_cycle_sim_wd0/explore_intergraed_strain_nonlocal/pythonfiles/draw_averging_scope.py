import numpy as np  
import matplotlib.pyplot as plt

xscope = 600

l1 = 100
x1 = np.linspace(-xscope, xscope, 1000)
y1 = np.exp(-4*x1**2/(l1**2))

l2 = 200
x2 = np.linspace(-xscope, xscope, 1000)
y2 = np.exp(-4*x2**2/(l2**2))

l3 = 300
x3 = np.linspace(-xscope, xscope, 1000)
y3 = np.exp(-4*x3**2/(l3**2))

plt.figure(figsize=(8, 4))
plt.plot(x1, y1, label='l = 100 m', color='blue')
plt.plot(x2, y2, label='l = 200 m', color='green')
plt.plot(x3, y3, label='l = 300 m', color='red')
plt.legend()
plt.xlabel('x (m)')
plt.ylabel('exp(-4x^2/l^2)')
plt.title('Nonlocal Averaging Scope with Different Length Scales')
plt.savefig('averaging_scope.png', dpi=300)
plt.show()