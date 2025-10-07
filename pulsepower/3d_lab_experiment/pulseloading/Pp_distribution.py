#!/usr/bin/env python3
"""Plot Pp as a function of Dsensor for several _fitting_param_alpha values."""

import numpy as np
import matplotlib.pyplot as plt

CONVERT_EFFICIENCY = 1.0
EM = 0.03
BASE_FACTORS = [9000]  # Different base factor multipliers
ALPHA_LIST = [0.35]
EXPONENT_LIST = [1.0, 0.5, 0.25, 0.1]


def compute_pp(
    alpha,
    dsensor_mm,
    base_factor,
    convert_efficiency=CONVERT_EFFICIENCY,
    em=EM,
    exponent=1.0,
):
    """Return Pp for the given alpha and sensor distance (Dsensor in mm)."""
    return (
        (0.1 * 1e6)
        * base_factor
        / np.pow(dsensor_mm, exponent)
        * (convert_efficiency * em) ** alpha
    )


def main():
    dsensor_min_mm = 1.6
    dsensor_max_mm = np.sqrt(1.6**2 + 30.0**2)
    dsensor_values = np.linspace(dsensor_min_mm, dsensor_max_mm, 400)

    plt.figure(figsize=(10, 6))

    # Plot for each base factor with different line styles
    linestyles = ["-", "--", "-."]
    for i, base_factor in enumerate(BASE_FACTORS):
        for alpha in sorted(ALPHA_LIST):
            for exponent in EXPONENT_LIST:
                pp_values = compute_pp(
                    alpha, dsensor_values, base_factor, exponent=exponent
                )
                label = f"base={base_factor}, alpha={alpha}, exponent={exponent}"
                plt.plot(
                    dsensor_values, pp_values, label=label, linestyle=linestyles[i]
                )

    plt.xlabel("Dsensor (mm)")
    plt.ylabel("Pp")
    plt.title("Pp vs Dsensor for different base factors and alpha values")
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.legend(title="Parameters", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xlim([1.6, 30.14])
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
