"""Matplotlib helpers: one consistent style and small plotting utilities (Agg backend)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, TypeAlias

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.ticker import LogFormatterSciNotation
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.floating[Any]]

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e1e0d9"

#: Fixed colour (and marker) per method, in a colour-blind-safe order; never cycled.
STYLE: dict[str, tuple[str, str]] = {
    "fd2": ("#2a78d6", "o"),
    "fd4": ("#eb6834", "s"),
    "numerov": ("#1baf7a", "^"),
    "shooting": ("#6250d6", "v"),
    "fourier": ("#e87ba4", "D"),
    "chebyshev": ("#e34948", "P"),
    "crank-nicolson": ("#2a78d6", "o"),
    "split-step": ("#eb6834", "s"),
    "rk4": ("#1baf7a", "^"),
    "direct": ("#898781", "x"),
    "cg": ("#2a78d6", "o"),
    "cg-jacobi": ("#6250d6", "v"),
    "cg-ic0": ("#eb6834", "s"),
    "mg": ("#1baf7a", "^"),
    "mg-cg": ("#e34948", "P"),
    "backward-euler": ("#eb6834", "s"),
}


def style(name: str) -> dict[str, object]:
    color, marker = STYLE.get(name, ("#52514e", "o"))
    return {"color": color, "marker": marker, "markersize": 5, "linewidth": 1.6}


def apply_theme() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": "#c3c2b7",
            "axes.labelcolor": MUTED,
            "axes.titlecolor": INK,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "text.color": INK,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "font.size": 9.5,
            "figure.dpi": 100,
        }
    )


def new_figure(ncols: int = 1, nrows: int = 1, width: float = 5.2, height: float = 4.0):
    apply_theme()
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows))
    return fig, axes


def loglog_slope_guide(
    ax, x: FloatArray, y0: float, slope: float, label: str, anchor: str = "min"
) -> None:
    """Dotted guide ``y = y0 (x / x_ref)^slope`` anchored at the min or max of ``x``."""
    xs = np.array([min(x), max(x)], dtype=float)
    x_ref = xs[1] if anchor == "max" else xs[0]
    ys = y0 * (xs / x_ref) ** slope
    ax.loglog(xs, ys, ls=":", color="#898781", lw=1.0)
    ax.annotate(
        label,
        (xs[-1], ys[-1]),
        textcoords="offset points",
        xytext=(-4, 6),
        ha="right",
        fontsize=8,
        color=MUTED,
    )


def save(fig: Figure, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    for ax in fig.axes:
        if ax.get_xscale() == "log":
            ax.xaxis.set_minor_formatter(LogFormatterSciNotation(minor_thresholds=(1.2, 0.3)))
        if ax.get_yscale() == "log":
            ax.yaxis.set_minor_formatter(LogFormatterSciNotation(minor_thresholds=(1.2, 0.3)))
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
