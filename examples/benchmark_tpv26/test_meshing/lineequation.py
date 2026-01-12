import numpy as np
import matplotlib.pyplot as plt

def compute_third_points_and_line_inequalities(x1, y1, x2, y2, angle_deg=15):
    p1 = np.array([x1, y1])
    p2 = np.array([x2, y2])
    # Midpoint and apex calculations
    m = (p1 + p2) / 2
    v = p2 - p1
    length = np.linalg.norm(v)
    v_unit = v / length
    v_perp = np.array([-v_unit[1], v_unit[0]])
    angle_rad = np.radians(angle_deg)
    h = (length / 2) * np.tan(angle_rad)
    p3a = m + h * v_perp
    p3b = m - h * v_perp

    def line_ineq(p_start, p_end, interior_pt):
        # Compute A, B, C for line through p_start->p_end
        x1_, y1_ = p_start
        x2_, y2_ = p_end
        A = y2_ - y1_
        B = x1_ - x2_
        C = x2_*y1_ - x1_*y2_
        val = A*interior_pt[0] + B*interior_pt[1] + C
        sign = ">" if val > 0 else "<"
        return A, B, C, sign

    solutions = {}
    for name, p3 in [("sol1", p3a), ("sol2", p3b)]:
        # Compute centroid as interior test point
        centroid = (p1 + p2 + p3) / 3
        ineqs = {}
        for side, (ps, pe) in [("P1-P2", (p1, p2)), ("P2-P3", (p2, p3)), ("P3-P1", (p3, p1))]:
            A, B, C, sign = line_ineq(ps, pe, centroid)
            # Format as "<coef>* x ± <coef>* y ± <const> sign 0"
            def term(coeff, var):
                return f"{abs(coeff):.6f} * {var}"
            parts = []
            parts.append(f"{'- ' if A < 0 else ''}{term(A, 'x')}")
            parts.append(f"{'- ' if B < 0 else '+ '}{term(B, 'y')}")
            parts.append(f"{'- ' if C < 0 else '+ '}{abs(C):.6f}")
            expr = " ".join(parts) + f" {sign} 0"
            ineqs[side] = expr
        solutions[name] = {"P3": p3, "inequalities": ineqs}

    return p1, p2, solutions

# Input points
x1, y1 = np.sqrt(3)/2 * 100 * 0.5, -0.5 * 100 * 0.5
x2, y2 = np.sqrt(3)/2 * 12000 * 1.03, -0.5 * 12000 * 1.03

# Compute
p1, p2, sols = compute_third_points_and_line_inequalities(x1, y1, x2, y2)

# Print results
print(f"P1: ({p1[0]:.6f}, {p1[1]:.6f})")
print(f"P2: ({p2[0]:.6f}, {p2[1]:.6f})\n")
for name, sol in sols.items():
    p3 = sol["P3"]
    print(f"{name.upper()} - P3: ({p3[0]:.6f}, {p3[1]:.6f})")
    for side, expr in sol["inequalities"].items():
        print(f"  {side}: {expr}")
    print()

# Plotting
plt.figure(figsize=(8,5))
plt.plot([p1[0], p2[0]], [p1[1], p2[1]], '-o', label='P1-P2')
for name, sol in sols.items():
    p3 = sol["P3"]
    plt.plot([p2[0], p3[0]], [p2[1], p3[1]], '--', label=f'P2-P3 {name}')
    plt.plot([p3[0], p1[0]], [p3[1], p1[1]], '--', label=f'P3-P1 {name}')
    plt.scatter(p3[0], p3[1], label=f'P3 {name}', s=80)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Triangle Solutions with Interior Inequalities')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()