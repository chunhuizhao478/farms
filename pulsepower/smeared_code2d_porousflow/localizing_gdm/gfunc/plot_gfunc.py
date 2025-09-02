import numpy as np
import matplotlib.pyplot as plt
import argparse

# g(d) = ((1−R) e^(−η d) + R − e^(−η)) / (1 − e^(−η))
def g_func(d, eta, R):
    num = (1.0 - R) * np.exp(-eta * d) + R - np.exp(-eta)
    den = 1.0 - np.exp(-eta)
    return num / den

def parse_float_list(s):
    return [float(x) for x in s.split(',')]

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Plot interaction function g(d) and l_eff/l.")
    p.add_argument("--eta", type=str, default="5", help="eta or comma list, e.g. '3,5,8'")
    p.add_argument("--R", type=str, default="0.005", help="R or comma list, e.g. '0.005,0.02'")
    p.add_argument("--l", type=float, default=1e-4, help="reference length l (only used for label)")
    p.add_argument("--out", type=str, default="g_plot.png", help="output figure file")
    args = p.parse_args()

    etas = parse_float_list(args.eta)
    Rs = parse_float_list(args.R)

    d = np.linspace(0.0, 1.0, 500)

    plt.figure(figsize=(9, 5))
    for eta in etas:
        for R in Rs:
            g = g_func(d, eta, R)
            leff_ratio = np.sqrt(np.clip(g, 0.0, None))  # l_eff / l
            plt.plot(d, g, lw=2, label=f"g(d), eta={eta}, R={R}")
            # plt.plot(d, leff_ratio, lw=1.5, ls="--", label=f"sqrt(g), eta={eta}, R={R}")

    plt.xlabel("damage d")
    plt.ylabel("value")
    plt.title(f"Interaction function g(d) and l_eff/l (l={args.l:g})")
    plt.grid(True, ls=":")
    plt.xlim(0, 1)
    plt.ylim(0, None)
    plt.legend(ncol=2)
    plt.tight_layout()
    plt.savefig(args.out, dpi=150)
    plt.show()