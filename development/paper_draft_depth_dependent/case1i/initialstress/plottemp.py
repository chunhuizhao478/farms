import numpy as np 
import matplotlib.pyplot as plt

a1_c = 24.335
a2_c = -0.11613

def T_constrained(z_val):
    return 279.1024 + 0.03982143 * (-z_val) - 0.000000559524 * (-z_val) * (-z_val)

z_examples = np.array([0e3, -5e3, -10e3, -15e3, -20e3], dtype=float)
T_examples = [T_constrained(zv) for zv in z_examples]
list(zip(z_examples, np.round(T_examples, 1)))

plt.figure()
plt.plot(z_examples,T_examples)
plt.show()