import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_path = "elasticity_csv.csv"

df = pd.read_csv(csv_path, comment='/', skip_blank_lines=True)

time_col = "time"
cols = [c for c in df.columns if c != time_col]

# Keep only data at or before 0.00006 s
t_max = 1e-5
df = df[df[time_col] <= t_max].sort_values(time_col)

t = df[time_col].to_numpy()

plt.figure(figsize=(10, 6))
for col in cols:
    y = df[col].to_numpy()
    m = np.isfinite(t) & np.isfinite(y)
    plt.plot(t[m], y[m], lw=1.0, label=col)

plt.xlabel("time (s)")
plt.ylabel("value")
plt.title("Energy Budget Time History")
plt.grid(True, ls=":", alpha=0.6)
plt.legend(loc="best", ncol=2)
plt.xlim(0, t_max)
# plt.yscale("symlog", linthresh=1e-6)

plt.tight_layout()
plt.show()