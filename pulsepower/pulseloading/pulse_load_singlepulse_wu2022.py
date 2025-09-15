import numpy as np
import matplotlib.pyplot as plt

"""Solve coupled channel pressure P(t) and volume V(t) growth for a single
plasma pulse using Wu et al. (2022) model.

Governing relations (notation follows user-provided snapshot):
  dWe/dt is prescribed (power deposition)
  dWe/dt = P dV/dt + 1/(gamma-1) d(PV)/dt
Re-arranged energy equation (used for dP/dt):
  (gamma-1) dWe/dt = gamma P dV/dt + V dP/dt  =>  dP/dt = ((gamma-1) dWe/dt - gamma P dV/dt)/V

Channel expansion law:
  dV/dt = 4/(n-1) * (pi * L * n * V)^{1/2} * psi^{1/(2n)} * rho0^{-1/2} * [ (P+psi)^{(n-1)/(2n)} - psi^{(n-1)/(2n)} ]

We assume cylindrical channel: V = pi * r_ch^2 * L  => r_ch = sqrt(V/(pi L))
"""

# ----------------------------- User / Physical Inputs -----------------------------
We = 30.0                 # J, total deposited electric energy
T = 2e-5                  # s, pulse duration for energy deposition
gamma = 1.33              # ratio of specific heats
n = 5                     # Wu et al. exponent for granite
rho0 = 1000               # kg/m^3, ambient density (only used for expansion model)
v0  = 1480                # m/s, speed of sound
psi = rho0 * v0 * v0 / n  # Pa, 
r0 = 1.6e-4               # m, initial channel radius guess (expansion model)
P0 = 0                    # Pa, initial channel pressure (start at ambient)

# ---------------------------- Derived / Helper Values ----------------------------
V0 = 4/3 * np.pi * r0**3
A = 3/2 * We / T  # W, peak power so that integral over pulse gives We (We = 2/3 A T)

def dWe_dt(t: float) -> float:
	"""Parabolic power deposition 0<=t<=T, zero otherwise.

	Peak chosen so integral = We (We = ∫ dWe/dt dt)
	Functional form: 4 A / T^2 * t (T - t)
	"""
	if 0.0 <= t <= T:
		return 4.0 * A / T**2 * t * (T - t)
	return 0.0

# Expansion law (spherical form):
# dV/dt = 2 (36π)^{1/3} sqrt(n)/(n-1) * rho0^{-1/2} * psi^{1/(2n)} * V^{2/3} * [ (P+psi)^{(n-1)/(2n)} - psi^{(n-1)/(2n)} ]
coeff_front = 2.0 * (36.0 * np.pi)**(1.0/3.0) * np.sqrt(n) / (n - 1.0) * rho0**(-0.5) * psi**(1.0/(2.0*n))
pow_exp = (n - 1.0) / (2.0 * n)

def dV_dt(P: float, V: float) -> float:
	"""Spherical channel volume growth rate (Wu et al. form provided in figure)."""
	V = max(V, 1e-30)
	PT_abs = max(P + psi, 1e-12)  # absolute pressure for exponent
	term = PT_abs**pow_exp - psi**pow_exp
	return coeff_front * (V**(2.0/3.0)) * term

def dP_dt(t: float, P: float, V: float) -> float:
	"""Pressure rate from rearranged energy equation."""
	Vc = max(V, 1e-30)
	dV = dV_dt(P, Vc)
	dwe = dWe_dt(t)
	return ((gamma - 1.0)*dwe - gamma * P * dV) / Vc

def rk4_step(t, P, V, dt):
	k1_V = dV_dt(P, V)
	k1_P = dP_dt(t, P, V)

	k2_V = dV_dt(P + 0.5*dt*k1_P, V + 0.5*dt*k1_V)
	k2_P = dP_dt(t + 0.5*dt, P + 0.5*dt*k1_P, V + 0.5*dt*k1_V)

	k3_V = dV_dt(P + 0.5*dt*k2_P, V + 0.5*dt*k2_V)
	k3_P = dP_dt(t + 0.5*dt, P + 0.5*dt*k2_P, V + 0.5*dt*k2_V)

	k4_V = dV_dt(P + dt*k3_P, V + dt*k3_V)
	k4_P = dP_dt(t + dt, P + dt*k3_P, V + dt*k3_V)

	V_next = V + dt*(k1_V + 2*k2_V + 2*k3_V + k4_V)/6.0
	P_next = P + dt*(k1_P + 2*k2_P + 2*k3_P + k4_P)/6.0
	return P_next, V_next

def simulate(t_end_factor: float = 8.0, n_steps: int = 4000):
	"""Full coupled P,V evolution (expanding radius)."""
	t_end = t_end_factor * T
	ts = np.linspace(0.0, t_end, n_steps)
	dt = ts[1] - ts[0]
	Ps = np.empty_like(ts)
	Vs = np.empty_like(ts)
	Ps[0] = P0
	Vs[0] = V0
	for i in range(n_steps - 1):
		Pn, Vn = rk4_step(ts[i], Ps[i], Vs[i], dt)
		Ps[i+1] = max(Pn, 0.0)
		Vs[i+1] = max(Vn, 1e-30)
	return ts, Ps, Vs

if __name__ == "__main__":
	# Coupled evolution only (fixed-radius code removed)
	ts, Ps, Vs = simulate()
	# Spherical radius from volume
	rchs = (3.0 * Vs / (4.0 * np.pi))**(1.0/3.0)

	# Power deposition
	plt.figure(figsize=(6,4))
	t_plot = np.linspace(0, ts[-1], 800)
	plt.plot(t_plot*1e6, [dWe_dt(t) for t in t_plot])
	plt.xlabel('Time ($\mu$s)')
	plt.ylabel('dWe/dt (W)')
	plt.title('Injected Power History')
	plt.grid(True)

	# Pressure vs time
	plt.figure(figsize=(6,4))
	plt.plot(ts*1e6, Ps*1e-6, label='Pressure')
	plt.xlabel('Time ($\mu$s)')
	plt.ylabel('P (MPa)')
	plt.title('Channel Pressure vs Time')
	plt.grid(True)
	plt.legend()

	# Radius vs time
	plt.figure(figsize=(6,4))
	plt.plot(ts*1e6, rchs*1e3, label='r_ch (spherical)')
	plt.xlabel('Time ($\mu$s)')
	plt.ylabel('$r_ch$ (mm)')
	plt.title('Channel Radius vs Time')
	plt.grid(True)
	plt.legend()

	# Pressure vs radius with threshold shading
	plt.figure(figsize=(6,4))
	r_vals_mm = rchs*1e3
	P_vals_MPa = Ps*1e-6
	r_threshold_mm = 1.6
	plt.plot(r_vals_mm, P_vals_MPa, label='P vs r_ch')
	cross_idx = np.where(r_vals_mm >= r_threshold_mm)[0]
	if cross_idx.size > 0:
		ic = cross_idx[0]
		if ic == 0:
			P_cross = P_vals_MPa[0]
		else:
			r1 = r_vals_mm[ic-1]; r2 = r_vals_mm[ic]
			P1 = P_vals_MPa[ic-1]; P2 = P_vals_MPa[ic]
			P_cross = P1 + (P2 - P1)*(r_threshold_mm - r1)/max(r2 - r1, 1e-30)
		plt.scatter([r_threshold_mm], [P_cross], color='red', zorder=5, label=f'r_ch={r_threshold_mm:.2f} mm')
		plt.axvline(r_threshold_mm, color='red', linestyle='--', alpha=0.8)
	if r_vals_mm.max() > r_threshold_mm:
		plt.axvspan(r_threshold_mm, r_vals_mm.max(), color='orange', alpha=0.18, label='r_ch ≥ 1.6 mm')
	plt.xlabel('$r_ch$ (mm)')
	plt.ylabel('P (MPa)')
	plt.title('Pressure vs Channel Radius')
	plt.grid(True)
	plt.legend()
	plt.tight_layout()

	# 2x2 combined summary figure
	fig, axes = plt.subplots(2, 2, figsize=(10,8))

	# (1) Power vs time
	ax = axes[0,0]
	ax.plot(t_plot*1e6, [dWe_dt(t) for t in t_plot])
	ax.set_xlabel('Time ($\\mu$s)')
	ax.set_ylabel('dWe/dt (W)')
	ax.set_title('Injected Power')
	ax.grid(True)

	# (2) Pressure vs time
	ax = axes[0,1]
	ax.plot(ts*1e6, Ps*1e-6, color = 'tab:blue')
	ax.set_xlabel('Time ($\\mu$s)')
	ax.set_ylabel('P (MPa)')
	ax.set_title('Pressure vs Time')
	ax.grid(True)

	# (3) Radius vs time
	ax = axes[1,0]
	ax.plot(ts*1e6, rchs*1e3, color = 'tab:blue')
	ax.set_xlabel('Time ($\\mu$s)')
	ax.set_ylabel('$r_{ch}$ (mm)')
	ax.set_title('Radius vs Time')
	ax.grid(True)

	# (4) Pressure vs radius with shading
	ax = axes[1,1]
	ax.plot(r_vals_mm, P_vals_MPa, label='P vs r_ch', color = 'tab:blue')
	ax.set_xlabel('$r_{ch}$ (mm)')
	ax.set_ylabel('P (MPa)')
	ax.set_title('Pressure vs Radius')
	ax.grid(True)
	if r_vals_mm.max() > r_threshold_mm:
		ax.axvspan(r_threshold_mm, r_vals_mm.max(), color='orange', alpha=0.18)
	cross_idx = np.where(r_vals_mm >= r_threshold_mm)[0]
	if cross_idx.size > 0:
		ic = cross_idx[0]
		if ic == 0:
			P_cross = P_vals_MPa[0]
		else:
			r1 = r_vals_mm[ic-1]; r2 = r_vals_mm[ic]
			P1 = P_vals_MPa[ic-1]; P2 = P_vals_MPa[ic]
			P_cross = P1 + (P2 - P1)*(r_threshold_mm - r1)/max(r2 - r1, 1e-30)
		ax.scatter([r_threshold_mm], [P_cross], color='red', zorder=5, label=f'r_ch={r_threshold_mm:.2f} mm')
		ax.axvline(r_threshold_mm, color='red', linestyle='--', alpha=0.8)
	ax.legend(loc='best')

	fig.tight_layout()
	fig.savefig('pulse_summary.png', dpi=300)
	fig.savefig('pulse_summary.pdf')
	plt.show()

	print(f"Coupled mode peak pressure: {Ps.max():.3e} Pa at t={ts[Ps.argmax()]:.3e} s")
