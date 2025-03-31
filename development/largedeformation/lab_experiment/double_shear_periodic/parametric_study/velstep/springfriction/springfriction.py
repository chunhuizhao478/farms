import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def system_odes(t, y, k, eps_po, A, theta0, B, eps_dot):
    """
    Returns the time derivatives [dsigma/dt, dtheta/dt]
    for the given state variables (sigma, theta).
    """
    sigma, theta = y

    # Equation for dsigma/dt
    dsigma_dt = k * (
        eps_dot 
        - eps_po * sigma * (theta / theta0)**(-B / A)
    )

    # Equation for dtheta/dt
    dtheta_dt = 1.0 - (
        eps_po * sigma * (theta / theta0)**(-B / A) * theta
    )

    return [dsigma_dt, dtheta_dt]


# ------------------------
# 1. Define parameters
# ------------------------
k       = 1.81e8       # example value
eps_po  = 5e-12      # example value for dot{ε}_{po}
A       = 2e7       # example value
theta0  = 1e4     # example value
B       = 1.2e7       # example value
eps_dot = 2e-4       # imposed strain rate (constant in this example)

# ------------------------
# 2. Initial conditions
# ------------------------
sigma0 = 0.0        # initial stress
theta0_val = 1e5    # initial "theta" (not to be confused with theta0 constant)
y0 = [sigma0, theta0_val]

# ------------------------
# 3. Time span for solution
# ------------------------
t_start = 0.0
t_end   = 100000.0
t_eval  = np.linspace(t_start, t_end, 400)

# ------------------------
# 4. Solve ODE system
# ------------------------
sol = solve_ivp(
    fun=lambda t, y: system_odes(t, y, k, eps_po, A, theta0, B, eps_dot),
    t_span=[t_start, t_end],
    y0=y0,
    t_eval=t_eval,
    method='RK45'
)

# Extract the solutions
t_vals    = sol.t
sigma_sol = sol.y[0, :]
theta_sol = sol.y[1, :]

# ------------------------
# 5. Plot results
# ------------------------
plt.figure(figsize=(10,4))

# --- sigma(t) ---
plt.subplot(1, 2, 1)
plt.plot(t_vals, sigma_sol, 'b-', label=r'$\sigma(t)$')
plt.xlabel('Time')
plt.ylabel(r'$\sigma$')
plt.title('Stress vs. Time')
plt.grid(True)
plt.legend()

# --- theta(t) ---
plt.subplot(1, 2, 2)
plt.plot(t_vals, theta_sol, 'r-', label=r'$\theta(t)$')
plt.xlabel('Time')
plt.ylabel(r'$\theta$')
plt.title('Theta vs. Time')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()