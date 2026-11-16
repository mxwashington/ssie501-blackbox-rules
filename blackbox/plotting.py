"""Shared figure style. One fixed colour per area (never re-ordered)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AREA_COLOR = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#eda100"}
INK = "#2b2b2b"
MUTED = "#6b6b6b"
GRID = "#e6e6e3"


def style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 160, "font.size": 10,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "axes.titlesize": 11,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "xtick.color": MUTED, "ytick.color": MUTED, "legend.frameon": False,
        "lines.linewidth": 2,
    })


style()
