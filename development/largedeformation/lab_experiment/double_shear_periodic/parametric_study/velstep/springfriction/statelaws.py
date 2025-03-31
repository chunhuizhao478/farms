import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Core parameters
C_g = 0.6        # Friction coefficient
tau = 1.0        # Normalized shear stress
N_ij = 1.0       # Normal stress
theta_o = 4e3    # Reference state
t_max = 100.0    # Total simulation time

# ===== New approach: Velocity-dependent constitutive parameters =====
# Define velocity-dependent A and B parameters
def A_v(v):
    """Velocity-dependent direct effect parameter"""
    A_low = 0.010  # A at low velocity (velocity weakening)
    A_high = 0.015  # A at high velocity (velocity strengthening)
    V_transition = 2.0  # Transition velocity (μm/s)
    width = 0.5  # Width of transition zone
    return A_low + (A_high - A_low) * (1 / (1 + np.exp(-(np.log10(v) - np.log10(V_transition))/width)))

def B_v(v):
    """Velocity-dependent evolution effect parameter"""
    B_low = 0.015  # B at low velocity (velocity weakening: B > A)
    B_high = 0.010  # B at high velocity (velocity strengthening: B < A)
    V_transition = 2.0  # Transition velocity (μm/s)
    width = 0.5  # Width of transition zone
    return B_low + (B_high - B_low) * (1 / (1 + np.exp(-(np.log10(v) - np.log10(V_transition))/width)))

# Define the plastic deformation rate function (Dp = V for simplicity)
def calc_Dp(V):
    return V

# Define state evolution laws with velocity-dependent parameters
def slip_law(t, theta, Dp):
    """Slip law for state evolution"""
    return -theta * Dp * np.log(theta * Dp / theta_o)

def steady_state_theta(Dp, A, B):
    """Calculate steady-state theta for given Dp, A, B"""
    return theta_o * (Dp)**(-B/A)

# Function to get steady-state friction vs. velocity curve
def get_steady_state_curve(velocities):
    mu_ss = []
    A_values = []
    B_values = []
    
    for v in velocities:
        A = A_v(v)
        B = B_v(v)
        A_values.append(A)
        B_values.append(B)
        
        theta_ss = steady_state_theta(v, A, B)
        mu = C_g * tau * N_ij * (theta_ss/theta_o)**(-B/A)
        mu_ss.append(mu)
        
    return mu_ss, A_values, B_values

# Function to simulate velocity step tests
def simulate_velocity_steps():
    # Define velocity step sequence (increasing steps)
    step_velocities = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
    step_times = np.linspace(10, 90, len(step_velocities))
    
    def velocity(t):
        """Return velocity at time t based on step sequence"""
        v = step_velocities[0]  # Default to first velocity
        for i, step_time in enumerate(step_times):
            if t >= step_time and (i == len(step_times)-1 or t < step_times[i+1]):
                v = step_velocities[i]
        return v
    
    # System of ODEs with velocity-dependent parameters
    def system(t, y):
        theta = y[0]
        
        V = velocity(t)
        Dp = calc_Dp(V)
        
        # Get velocity-dependent parameters
        A = A_v(V)
        B = B_v(V)
        
        # Calculate friction coefficient
        mu = C_g * tau * N_ij * (theta/theta_o)**(-B/A)
        
        # Use slip law for state evolution
        dtheta_dt = slip_law(t, theta, Dp)
            
        return [dtheta_dt, mu, A, B]  # Return state variable, friction, and parameters
    
    # Initial conditions: start with steady-state theta for first velocity
    V_init = step_velocities[0]
    A_init = A_v(V_init)
    B_init = B_v(V_init)
    theta_init = steady_state_theta(V_init, A_init, B_init)
    y0 = [theta_init, 0.0, A_init, B_init]  # [theta, mu, A, B]
    
    # Time points
    t_span = (0, t_max)
    t_eval = np.linspace(0, t_max, 2000)
    
    # Solve ODE system
    solution = solve_ivp(system, t_span, y0, method='RK45', t_eval=t_eval)
    
    return solution.t, solution.y, [velocity(t) for t in solution.t], step_times, step_velocities

# Function to simulate continuous velocity ramp
def simulate_velocity_ramp():
    # Define velocity ramp
    V_min = 0.1
    V_max = 10.0
    
    def velocity(t):
        """Return velocity at time t based on logarithmic ramp"""
        if t < 10:
            return V_min
        elif t > 90:
            return V_max
        else:
            # Log-linear ramp from V_min to V_max
            t_normalized = (t - 10) / 80  # normalize to 0-1
            log_v = np.log10(V_min) + t_normalized * (np.log10(V_max) - np.log10(V_min))
            return 10**log_v
    
    # System of ODEs with velocity-dependent parameters
    def system(t, y):
        theta = y[0]
        
        V = velocity(t)
        Dp = calc_Dp(V)
        
        # Get velocity-dependent parameters
        A = A_v(V)
        B = B_v(V)
        
        # Calculate friction coefficient
        mu = C_g * tau * N_ij * (theta/theta_o)**(-B/A)
        
        # Use slip law for state evolution
        dtheta_dt = slip_law(t, theta, Dp)
            
        return [dtheta_dt, mu, A, B]  # Return state variable, friction, and parameters
    
    # Initial conditions: start with steady-state theta for first velocity
    V_init = V_min
    A_init = A_v(V_init)
    B_init = B_v(V_init)
    theta_init = steady_state_theta(V_init, A_init, B_init)
    y0 = [theta_init, 0.0, A_init, B_init]  # [theta, mu, A, B]
    
    # Time points
    t_span = (0, t_max)
    t_eval = np.linspace(0, t_max, 2000)
    
    # Solve ODE system
    solution = solve_ivp(system, t_span, y0, method='RK45', t_eval=t_eval)
    
    return solution.t, solution.y, [velocity(t) for t in solution.t]

# Run the simulations
print("Running simulations...")

# Calculate steady-state friction vs. velocity curve
V_range = np.logspace(-1, 1, 100)  # Range of velocities from 0.1 to 10 μm/s
mu_ss, A_values, B_values = get_steady_state_curve(V_range)

# Run velocity step simulation
t_step, y_step, v_step, step_times, step_velocities = simulate_velocity_steps()
theta_step = y_step[0]
mu_step = y_step[1]
A_step = y_step[2]
B_step = y_step[3]

# Run velocity ramp simulation
t_ramp, y_ramp, v_ramp = simulate_velocity_ramp()
theta_ramp = y_ramp[0]
mu_ramp = y_ramp[1]
A_ramp = y_ramp[2]
B_ramp = y_ramp[3]

# ====== Create plots ======

# Figure 1: Steady-state friction vs. velocity curve
fig1, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Plot steady-state friction
ax1.semilogx(V_range, mu_ss, 'k-', linewidth=2)
ax1.set_ylabel('Steady-State Friction Coefficient')
ax1.set_title('Steady-State Friction vs. Velocity')
ax1.grid(True)

# Add marks for velocity strengthening/weakening regions
transition_v = 2.0
ax1.axvline(x=transition_v, color='r', linestyle='--', 
            label=f'Transition Velocity ({transition_v} μm/s)')
ax1.text(0.2, 0.95*max(mu_ss), 'Velocity Weakening', 
         horizontalalignment='center', color='darkred')
ax1.text(7.0, 0.95*max(mu_ss), 'Velocity Strengthening', 
         horizontalalignment='center', color='darkblue')
ax1.legend()

# Plot A and B values
ax2.semilogx(V_range, A_values, 'b-', linewidth=2, label='A(v)')
ax2.semilogx(V_range, B_values, 'g-', linewidth=2, label='B(v)')
ax2.set_xlabel('Velocity (μm/s)')
ax2.set_ylabel('Parameter Value')
ax2.set_title('Velocity-Dependent Parameters')
ax2.grid(True)
ax2.legend()

plt.tight_layout()
plt.savefig('steady_state_curve_with_parameters.png', dpi=300)

# Figure 2: Velocity step test results
fig2, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

# Plot applied velocity
ax1.step(t_step, v_step, 'k-', where='post', linewidth=2)
ax1.set_ylabel('Velocity (μm/s)')
ax1.set_title('Velocity Step Test')
ax1.grid(True)
for time in step_times:
    ax1.axvline(x=time, color='r', linestyle='--', alpha=0.3)

# Plot friction coefficient
ax2.plot(t_step, mu_step, 'b-', linewidth=2)
ax2.set_ylabel('Friction Coefficient')
ax2.grid(True)

# Plot state variable
ax3.plot(t_step, theta_step, 'g-', linewidth=2)
ax3.set_ylabel('State Variable θ')
ax3.set_yscale('log')
ax3.grid(True)

# Plot A/B ratio
ab_ratio = A_step / B_step
ax4.plot(t_step, ab_ratio, 'r-', linewidth=2)
ax4.axhline(y=1.0, color='k', linestyle='--', 
            label='A/B = 1 (Transition Line)')
ax4.set_ylabel('A/B Ratio')
ax4.set_xlabel('Time (s)')
ax4.grid(True)
ax4.legend()

plt.tight_layout()
plt.savefig('velocity_step_test.png', dpi=300)

# Figure 3: Velocity ramp results
fig3, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

# Plot applied velocity
ax1.plot(t_ramp, v_ramp, 'k-', linewidth=2)
ax1.set_ylabel('Velocity (μm/s)')
ax1.set_title('Velocity Ramp Test')
ax1.set_yscale('log')
ax1.grid(True)

# Plot friction coefficient
ax2.plot(t_ramp, mu_ramp, 'b-', linewidth=2)
ax2.set_ylabel('Friction Coefficient')
ax2.grid(True)

# Plot state variable
ax3.plot(t_ramp, theta_ramp, 'g-', linewidth=2)
ax3.set_ylabel('State Variable θ')
ax3.set_yscale('log')
ax3.grid(True)

# Plot A/B ratio
ab_ratio_ramp = A_ramp / B_ramp
ax4.plot(t_ramp, ab_ratio_ramp, 'r-', linewidth=2)
ax4.axhline(y=1.0, color='k', linestyle='--', 
            label='A/B = 1 (Transition Line)')
ax4.set_ylabel('A/B Ratio')
ax4.set_xlabel('Time (s)')
ax4.grid(True)
ax4.legend()

plt.tight_layout()
plt.savefig('velocity_ramp_test.png', dpi=300)

# Figure 4: Friction vs. Velocity during ramp (with steady-state curve comparison)
fig4, ax = plt.subplots(figsize=(10, 8))

# Plot steady-state curve
ax.semilogx(V_range, mu_ss, 'k-', linewidth=2, label='Steady-State Curve')

# Plot actual friction vs. velocity during ramp (may show hysteresis)
ax.semilogx(v_ramp, mu_ramp, 'r--', linewidth=1.5, alpha=0.7, label='Simulation Result')

# Mark transition from weakening to strengthening
ax.axvline(x=2.0, color='b', linestyle='--', 
           label='Transition Velocity (2.0 μm/s)')

ax.set_xlabel('Velocity (μm/s)')
ax.set_ylabel('Friction Coefficient')
ax.set_title('Friction-Velocity Relationship')
ax.grid(True)
ax.legend()

plt.tight_layout()
plt.savefig('friction_velocity_comparison.png', dpi=300)

print("Simulation complete. Key observations:")

# Calculate if each velocity step shows weakening or strengthening
for i in range(len(step_velocities)-1):
    v1 = step_velocities[i]
    v2 = step_velocities[i+1]
    
    # Find indices right before and after step
    t1 = step_times[i] - 0.1
    t2 = step_times[i] + 5  # Allow time to reach near steady state
    
    idx1 = np.argmin(np.abs(t_step - t1))
    idx2 = np.argmin(np.abs(t_step - t2))
    
    mu1 = mu_step[idx1]
    mu2 = mu_step[idx2]
    
    behavior = "strengthening" if mu2 > mu1 else "weakening"
    print(f"Step from {v1} to {v2} μm/s shows velocity {behavior} (Δμ = {mu2-mu1:.5f})")

# Report A/B ratios at different velocities
print("\nA/B ratios at different velocities:")
test_velocities = [0.2, 1.0, 2.0, 5.0, 8.0]
for v in test_velocities:
    A = A_v(v)
    B = B_v(v)
    expected = "strengthening" if A > B else "weakening"
    print(f"At v = {v} μm/s: A/B = {A/B:.3f} (expects velocity {expected})")