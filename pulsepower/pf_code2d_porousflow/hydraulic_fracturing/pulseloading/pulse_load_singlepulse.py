import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


def pressure(t, alpha, beta, t0, p0=1.0):
    """Normalized pressure history (p(t0)=p0, p→0 as t→t0+td)."""
    # Vectorized operations; assume alpha<beta
    num = np.exp(-alpha * t) - np.exp(-beta * t)
    den = np.exp(-alpha * t0) - np.exp(-beta * t0)
    return p0 * num / den


def peak_time(alpha, beta):
    """Return analytical peak time where dp/dt=0 for the bi-exponential difference."""
    return np.log(beta / alpha) / (beta - alpha)


def objective(params, t0, td, epsilon):
    """Scalar objective enforcing peak at t0 and near-zero at t0+td."""
    alpha, beta = params
    if not (alpha > 0 and beta > 0 and alpha < beta):
        return 1e9
    pt_err = abs(t0 - peak_time(alpha, beta))
    decay_val = pressure(t0 + td, alpha, beta, t0)
    decay_err = abs(decay_val - epsilon)
    # Weighted sum (can tune weights if needed)
    return pt_err + decay_err


def solve_params(t0, td, epsilon):
    # Heuristic initial guesses
    alpha_guess = 1.0 / td
    beta_guess = max(5.0 / t0, alpha_guess * 1.1)
    x0 = [alpha_guess, beta_guess]
    res = minimize(
        objective,
        x0,
        args=(t0, td, epsilon),
        method="Powell",
        options={"xatol": 1e-12, "fatol": 1e-12, "maxiter": 20000, "disp": False},
    )
    if not res.success:
        print(
            f"WARNING: optimization did not fully converge for t0={t0:.2e}s: {res.message}"
        )
    alpha_opt, beta_opt = res.x
    return alpha_opt, beta_opt, res


def main():
    # User settings (time unit: seconds)
    td = 0.0025  # fixed decay time (0.1 s)
    t0_list = [5e-5]  # peak time (5 ms = 0.005 s)
    epsilon = 1e-6  # target near-zero value at t0+td
    colors = ["tab:blue", "tab:red"]

    # Time range covers entire decay for the largest t0
    t_end = max(t0_list) + td
    t_values = np.linspace(0.0, t_end, 1500)

    # Styling aligned with wu2022 script
    plt.figure(figsize=(6, 4))
    for t0, color in zip(t0_list, colors):
        alpha_opt, beta_opt, res = solve_params(t0, td, epsilon)
        p_vals = pressure(t_values, alpha_opt, beta_opt, t0)
        label = rf"$t_0={t0:.3f}\,s$"
        plt.plot(t_values, p_vals, color=color, lw=2, label=label)
        # Mark peak
        plt.axvline(t0, color=color, ls="--", lw=1, alpha=0.6)
        print(
            f"t0={t0:.3e} s: alpha={alpha_opt:.3e}, beta={beta_opt:.3e}, peak_time_err={abs(t0 - peak_time(alpha_opt, beta_opt)):.2e}, p(t0+td)={pressure(t0 + td, alpha_opt, beta_opt, t0):.2e}"
        )

    plt.xlabel("Time (s)")
    plt.ylabel("$p(t)/p_{0}$")
    plt.title("Normalized Pressure Profiles")
    plt.grid(True, which="both", alpha=0.35)
    plt.xlim(0, t_end)
    plt.ylim(0, 1.05)
    plt.legend(frameon=True, fontsize=9)
    plt.tight_layout()
    plt.savefig("singlepulse_profiles.png", dpi=300)
    plt.savefig("singlepulse_profiles.pdf")
    plt.show()


if __name__ == "__main__":
    main()
