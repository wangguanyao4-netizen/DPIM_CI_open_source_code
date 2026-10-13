"""Plot-only functions extracted from the accepted internal-n renderer."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIG_DIR = OUT_DIR = Path("outputs")
METHODS = ["t distribution", "percentile bootstrap", "bootstrap-t"]


METHOD_LABEL = {
    "t distribution": "Student t",
    "percentile bootstrap": "percentile bootstrap",
    "bootstrap-t": "bootstrap-t",
}


METHOD_COLOR = {
    "t distribution": "#245DA0",
    "percentile bootstrap": "#C25746",
    "bootstrap-t": "#2F7D45",
}


N_VALUES = [192, 384, 768]


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.edgecolor": "#20252B",
            "axes.linewidth": 0.75,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": False,
            "ytick.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def plot_n_coverage(df: pd.DataFrame, model: str, output_name: str) -> None:
    frame = df[df["model"] == model].copy()
    fig, axes = plt.subplots(len(N_VALUES), len(METHODS), figsize=(8.9, 5.65), sharex=True, sharey=True)
    for i, n in enumerate(N_VALUES):
        for j, method in enumerate(METHODS):
            ax = axes[i, j]
            sub = frame[(frame["n"] == n) & (frame["method"] == method)].sort_values("y")
            color = METHOD_COLOR[method]
            ax.plot(sub["y"], sub["coverage"], color=color, linewidth=1.05)
            selected = np.unique(np.round(np.linspace(0, len(sub) - 1, 17)).astype(int))
            ax.scatter(
                sub["y"].to_numpy()[selected],
                sub["coverage"].to_numpy()[selected],
                s=13,
                facecolor="white",
                edgecolor=color,
                linewidth=0.75,
                zorder=3,
            )
            ax.axhline(0.95, color="#5C626B", linestyle=":", linewidth=0.8)
            if i == 0:
                ax.set_title(METHOD_LABEL[method], fontsize=10.0, pad=3)
            if j == 0:
                ax.set_ylabel(f"$n={n}$\ncoverage", fontsize=9.1)
            if i == len(N_VALUES) - 1:
                ax.set_xlabel(r"response coordinate $y$", fontsize=9.1)
            ax.tick_params(labelsize=8.4)
            ax.grid(True, axis="y", color="#D7DCE2", linewidth=0.35, alpha=0.70)
    lo=min(0.88,float(frame['coverage'].min())-.01)
    hi=max(0.985,float(frame['coverage'].max())+.01)
    axes[0, 0].set_ylim(max(0,lo),min(1,hi))
    fig.tight_layout(w_pad=0.85, h_pad=0.70)
    for ext in ("pdf",):
        fig.savefig(FIG_DIR / f"{output_name}.{ext}", dpi=320, bbox_inches="tight", pad_inches=0.035)
    plt.close(fig)

