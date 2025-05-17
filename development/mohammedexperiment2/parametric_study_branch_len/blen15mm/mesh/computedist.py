import numpy as np

# Given points (x, y, z)
P9  = np.array([0.141427, 0.116669, 0.0])
P10 = np.array([0.112338, 0.169146, 0.0])

# Compute the unit vector from P9 to P10
direction = P10 - P9
unit_dir  = direction / np.linalg.norm(direction)

# Distances from P9 along the P9→P10 line
distances = [0.01, 0.015, 0.06]

# Compute and print the new points
for d in distances:
    P_new = P9 + d * unit_dir
    print(f"Point at distance {d:.3f} from P9: ({P_new[0]:.6f}, {P_new[1]:.6f}, {P_new[2]:.6f})")
