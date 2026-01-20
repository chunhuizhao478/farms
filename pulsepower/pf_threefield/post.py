import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

csv_path = "elasticity_phasefield_csv.csv"

df = pd.read_csv(csv_path, comment="/", skip_blank_lines=True)

time_col = "time"
# Columns to hide (do not plot)
hidden_cols = {"damping_work"}


def format_label(col_name):
    """Convert underscore_separated_name to Title Case Name"""
    return " ".join(word.capitalize() for word in col_name.split("_"))


# Only plot columns that are not time and not in hidden list
cols = [c for c in df.columns if c != time_col and c not in hidden_cols]

# Keep only data at or before 0.00006 s
t_max = 2e-5
df = df[df[time_col] <= t_max].sort_values(time_col)

t = df[time_col].to_numpy() * 1e6  # Convert to microseconds

plt.figure(figsize=(10, 6))
for col in cols:
    y = df[col].to_numpy()
    m = np.isfinite(t) & np.isfinite(y)
    plt.plot(t[m], y[m], lw=1.0, label=format_label(col))

plt.xlabel("time ($\mu s$)", fontsize=18)
plt.ylabel("Energy (J)", fontsize=18)
plt.title("Energy Budget Time History", fontsize=20)
plt.grid(True, ls=":", alpha=0.6)
plt.legend(loc="best", ncol=2)
plt.xlim(0, t_max * 1e6)
# plt.ylim(0, 30)

plt.tight_layout()
plt.savefig("energy_budget.png", dpi=300, bbox_inches="tight")
plt.show()
