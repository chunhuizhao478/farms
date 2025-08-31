import matplotlib.pyplot as plt
from matplotlib import rcParams


def plot_tpv32_table(out_path: str = "tpv32_velocity_structure_table.png") -> None:
    """Render the TPV32 Velocity Structure table and save as an image.

    Parameters
    ----------
    out_path : str
        Output image filename (PNG recommended).
    """

    # Table content from TPV32 spec
    header = [
        "Depth = y (meters)",
        r"$V_P$ (meter/second)",
        r"$V_S$ (meter/second)",
        r"Density = $\rho$ (kg/m$^3$)",
    ]

    rows = [
        ["0", "2200", "1050", "2200"],
        ["500", "3000", "1400", "2450"],
        ["1000", "3600", "1950", "2550"],
        ["1600", "4400", "2500", "2600"],
        ["2400", "4800", "2800", "2600"],
        ["3600", "5250", "3100", "2620"],
        ["5000", "5500", "3250", "2650"],
        ["9000", "5750", "3450", "2720"],
        ["11000", "6100", "3600", "2750"],
        ["15000", "6300", "3700", "2900"],
        ["All depths y > 15000", "6300", "3700", "2900"],
    ]

    # Figure setup
    rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10,
    })
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.axis("off")

    # Title
    title = "TPV32 Velocity Structure"
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)

    # Build table
    table = ax.table(
        cellText=rows,
        colLabels=header,
        loc="center",
        cellLoc="center",
    )

    # Styling: column widths and cell fonts
    # Slightly wider first column for the final row text
    # Matplotlib indexes columns by (row, col) internally; use auto width adjustments
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.3)  # widen columns a bit, increase row height

    # Bold header
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#eaeaea")

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    plot_tpv32_table()
