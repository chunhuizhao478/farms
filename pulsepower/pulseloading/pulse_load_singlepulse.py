import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def pressure(t, alpha, beta, t0, p0=1.0):
    """Compute pressure at time t"""
    numerator = np.exp(-alpha * t) - np.exp(-beta * t)
    denominator = np.exp(-alpha * t0) - np.exp(-beta * t0)
    return p0 * numerator / denominator

def objective(params, t0, td, epsilon=1e-6):
    """
    Objective function to minimize:
    1. Error in peak time condition
    2. Error in decay condition
    """
    alpha, beta = params
    
    # Avoid invalid parameter combinations
    if alpha >= beta or alpha <= 0 or beta <= 0:
        return 1e6
    
    # Peak time condition error
    peak_time_error = np.abs(t0 - (1/(beta - alpha)) * np.log(beta/alpha))
    
    # Decay condition error
    p_at_decay = pressure(t0 + td, alpha, beta, t0)
    decay_error = np.abs(p_at_decay - epsilon)
    
    return peak_time_error + decay_error

# Parameters for both cases
td_values = [5e-5]  # Fixed decay time
t0 = 2e-6  # Two different rise times
epsilon = np.finfo(float).eps 

# Colors for the plots
colors = ['blue', 'red']
labels = ['td = 5e-5 s', 'td = 1e-4 s']

# Create plot
plt.figure(figsize=(10, 6))

# Generate curves for each t0
for td, color, label in zip(td_values, colors, labels):
    # Initial guess based on characteristic times
    alpha_guess = 1.0/td
    beta_guess = 5.0/t0
    initial_guess = [alpha_guess, beta_guess]
    
    # Optimization
    result = minimize(objective, initial_guess, 
                    args=(t0, td, epsilon),
                    method='powell',
                    options={'xatol': 1e-12, 'fatol': 1e-12, 'maxiter': 50000})
    
    if result.success:
        alpha_opt, beta_opt = result.x
        print(f"\nResults for {label}:")
        print(f"alpha = {alpha_opt:.3e}, beta = {beta_opt:.3e}")
        
        # Generate and plot pressure curve
        t_values = np.linspace(0, 1e-4, 1000)
        p_values = pressure(t_values, alpha_opt, beta_opt, t0)
        plt.plot(t_values, p_values, color=color, linewidth=2, label=label)
        
        # Print verification
        print("Verification:")
        print(f"Peak time error: {np.abs(t0 - (1/(beta_opt - alpha_opt)) * np.log(beta_opt/alpha_opt)):.2e}")
        print(f"Value at t0 + td: {pressure(t0 + td, alpha_opt, beta_opt, t0):.2e}")
    else:
        print(f"Optimization failed for {label}:", result.message)

plt.xlabel('Time (s)')
plt.ylabel('p(t)/p₀')
plt.title(f'Pressure Profiles (td = {td:.1e}s)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.xlim(0, 2e-5)
plt.ylim(0, 1.1)
plt.tight_layout()
plt.show()