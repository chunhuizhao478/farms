#!/usr/bin/env python3
"""Plot Pp as a function of Dsensor for several _fitting_param_alpha values."""

import numpy as np
import matplotlib.pyplot as plt

CONVERT_EFFICIENCY = 1.0
EM = 0.6
BASE_FACTOR = (0.1 * 1e6) * 9000
ALPHA_LIST = [0.01, 0.05, 0.1, 0.2, 0.35]


def compute_pp(alpha, dsensor_mm, convert_efficiency=CONVERT_EFFICIENCY, em=EM):
    """Return Pp for the given alpha and sensor distance (Dsensor in mm)."""
    return BASE_FACTOR / dsensor_mm * (convert_efficiency * em) ** alpha


def main():
    dsensor_min_mm = 1.6
    dsensor_max_mm = np.sqrt(1.6 ** 2 + 30.0 ** 2)
    dsensor_values = np.linspace(dsensor_min_mm, dsensor_max_mm, 400)

    plt.figure(figsize=(8, 5))

    for alpha in sorted(ALPHA_LIST):
        pp_values = compute_pp(alpha, dsensor_values)
        label = f"alpha = {alpha}"
        plt.plot(dsensor_values, pp_values, label=label)

    plt.xlabel("Dsensor (mm)")
    plt.ylabel("Pp")
    plt.title("Pp vs Dsensor for selected _fitting_param_alpha values")
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.legend(title="_fitting_param_alpha")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
