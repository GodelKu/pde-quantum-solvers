"""Run every study and write PNG figures, ``results.json`` and ``tables.md`` (``sbench report``).

Nothing here is hand-typed: the studies write their numbers to ``docs/results/results.json``,
and every README table is rendered from that JSON by :mod:`schrodinger_bench.tables` into the
blocks between ``<!-- BEGIN:name -->`` / ``<!-- END:name -->`` markers. Wall times are the median of
:data:`schrodinger_bench.timing.REPEATS` runs after one warm-up run, with their range.
"""

from __future__ import annotations

import json
import platform
import time
import warnings
from collections.abc import Callable
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.linalg import eig_banded, eigh_tridiagonal
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar

from . import plotting as pl
from .convergence import (
    HARD_WALL_SIZES,
    POTENTIALS,
    STATIONARY_SIZES,
    TIME_STEPS,
    ConvergenceStudy,
    conservation_comparison,
    rk4_stability_demo,
    spatial_study,
    stationary_study,
    temporal_study,
)
from .elliptic import (
    POISSON_SOLVERS,
    heat_exact_semidiscrete,
    heat_solve,
    l2_error,
    manufactured_problem,
    solve_poisson,
)
from .grid import Grid
from .potentials import DoubleWell, Harmonic, Morse, RectangularBarrier
from .stationary import SOLVERS, align_sign, fd4_banded_matrix
from .tables import fill, md_table, num, render_tables, sci, tables_markdown
from .timedep import (
    TIME_SOLVERS,
    CrankNicolson,
    SplitStepFourier,
    coherent_state_nd,
    gaussian_packet,
)
from .timing import REPEATS, repeat_runs, time_stats
from .tunneling import packet_transmission, packet_transmission_later, transmission_fd

__all__ = ["environment", "format_study", "generate_report", "md_table", "num", "sci"]


def environment() -> dict[str, str]:
    """Library versions and the OS family only (no host names, kernel strings or hardware)."""
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "os": platform.system(),
    }


def update_readme(readme: Path, tables: dict[str, str]) -> list[str]:
    """Replace the text between BEGIN/END markers; returns the names that were found."""
    text, found = fill(readme.read_text(encoding="utf-8"), tables)
    readme.write_text(text, encoding="utf-8")
    return found


def _clean(obj: Any) -> Any:
    """JSON-safe copy: numpy scalars to Python numbers, NaN to None, +-inf to "inf"/"-inf"."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items()}
    if isinstance(obj, list | tuple | np.ndarray):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, float | np.floating):
        if obj != obj:
            return None
        return float(obj) if np.isfinite(obj) else ("inf" if obj > 0 else "-inf")
    return obj


def _clock(fn: Callable[[], object], repeats: int) -> dict[str, float]:
    """Repeated wall time of a call that has no timer of its own."""

    def run() -> float:
        t0 = time.perf_counter()
        fn()
        return time.perf_counter() - t0

    return time_stats(repeat_runs(run, repeats))


def study_to_dict(s: ConvergenceStudy) -> dict[str, Any]:
    """Errors and orders of a study.

    Single-run timings are left out; see :mod:`schrodinger_bench.timing`.
    """
    rows = [{k: v for k, v in r.items() if k != "seconds"} for r in s.rows()]
    return {
        "name": s.name,
        "parameter": s.parameter,
        "reference": s.reference,
        "expected_order": s.expected,
        "fitted_order": s.fitted_order,
        "rows": rows,
    }


def format_study(s: ConvergenceStudy) -> str:
    """Plain-text table for the terminal."""
    head = f"{s.name}   (error vs {s.reference}; nominal order: {s.expected})"
    cols = [s.parameter, "size", "error", "observed order"]
    lines = [head, "  " + "  ".join(f"{c:>14s}" for c in cols)]
    orders = [float("nan"), *s.pairwise_orders.tolist()]
    for i in range(len(s.errors)):
        size = str(s.sizes[i]) if s.sizes else ""
        lines.append(
            "  "
            + "  ".join(
                f"{v:>14s}" for v in (sci(s.steps[i], 3), size, sci(s.errors[i], 3), num(orders[i]))
            )
        )
    if not s.is_spectral:
        lines.append(f"  least-squares fitted order: {s.fitted_order:.3f}")
    else:
        lines.append("  spectral method: error decays faster than any power of h")
    return "\n".join(lines)


# --------------------------------------------------------------------------------------------
# 1. Stationary problem
# --------------------------------------------------------------------------------------------


def stationary_section(fig_dir: Path, quick: bool, repeats: int) -> dict[str, Any]:
    studies: dict[str, Any] = {}
    fig, axes = pl.new_figure(ncols=3, width=5.0, height=4.0)
    for ax, (pname, pot) in zip(axes, POTENTIALS.items(), strict=True):
        studies[pname] = {}
        for method, solver in SOLVERS.items():
            sizes = STATIONARY_SIZES[method][:4] if quick else STATIONARY_SIZES[method]
            st = stationary_study(method, pname, sizes=sizes)
            entry = study_to_dict(st)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                runs = repeat_runs(
                    partial(solver.solve, pot, pot.box[0], pot.box[1], sizes[-1], 4), repeats
                )
            entry["time_largest_n"] = time_stats([r.wall_time for r in runs])
            studies[pname][method] = entry
            err = np.maximum(st.errors, 1e-16)
            ax.loglog(np.array(st.sizes), err, label=method, **pl.style(method))
        ax.set_title(f"{pname}  (reference: {studies[pname]['fd2']['reference']})")
        ax.set_xlabel("N (interior nodes / modes)")
        ax.set_ylabel("max error of lowest 4 levels")
        ax.set_ylim(1e-15, 1e1)
        if pname == "harmonic":
            n_ref = np.array([40.0, 1280.0])
            pl.loglog_slope_guide(ax, n_ref, 3e-1, -2, "N^-2")
            pl.loglog_slope_guide(ax, n_ref, 5e-3, -4, "N^-4")
        ax.legend(loc="lower left", ncol=2)
    pl.save(fig, fig_dir / "stationary_convergence.png")
    data: dict[str, Any] = {"studies": studies}

    # Wavefunction accuracy against the exact oscillator eigenfunctions, plus eigenfunction plot
    ho = Harmonic()
    wf: dict[str, Any] = {}
    for method, solver in SOLVERS.items():
        res = solver.solve(ho, *ho.box, 200, 6)
        errs = []
        for k in (0, 3):
            ex = ho.exact_state(res.x, k)
            psi = align_sign(res.states[:, k], ex)
            errs.append(float(np.max(np.abs(psi - ex))))
        wf[method] = {
            "n": 200,
            "psi0": errs[0],
            "psi3": errs[1],
            "energy_6_levels": float(np.max(np.abs(res.energies - ho.exact_energies(6)))),
        }
    data["wavefunctions"] = wf

    fig, axes = pl.new_figure(ncols=3, width=5.0, height=4.2)
    for ax, (pname, pot) in zip(axes, POTENTIALS.items(), strict=True):
        x_min, x_max = pot.box
        res = SOLVERS["chebyshev"].solve(pot, x_min, x_max, 200, 4)
        xs = np.linspace(x_min, x_max, 600)
        ax.plot(xs, pot(xs), color=pl.INK, lw=1.4)
        scale = 0.45 * max(float(np.diff(res.energies).min()), 1.0)
        for k in range(4):
            c, _ = pl.STYLE[list(pl.STYLE)[k]]
            ax.axhline(res.energies[k], color=c, lw=0.6, alpha=0.6)
            sign = 1.0 if res.states[:, k][np.argmax(np.abs(res.states[:, k]))] > 0 else -1.0
            ax.plot(
                res.x,
                res.energies[k] + sign * scale * res.states[:, k] / np.max(np.abs(res.states)),
                color=c,
                lw=1.5,
                label=f"n={k}",
            )
        top = res.energies[-1] + 3.2
        ax.set_ylim(min(float(pot(xs).min()), 0.0) - 0.5, top)
        lim = {"harmonic": 5.0, "double-well": 4.5, "morse": None}[pname]
        ax.set_xlim(-lim, lim) if lim else ax.set_xlim(-1.5, 8.0)
        ax.set_title(pname)
        ax.set_xlabel("x")
        ax.set_ylabel("E")
        ax.legend(loc="upper center", ncol=4)
    pl.save(fig, fig_dir / "stationary_states.png")

    # Double-well tunnel splitting from several methods
    dw = DoubleWell()
    split: dict[str, Any] = {}
    for method in ("fd2", "fd4", "numerov", "shooting", "fourier", "chebyshev"):
        n = 400 if method in ("fd2", "fd4", "numerov", "shooting") else 96
        r = SOLVERS[method].solve(dw, *dw.box, n, 2)
        split[method] = {
            "n": n,
            "e0": float(r.energies[0]),
            "splitting": float(r.energies[1] - r.energies[0]),
        }
    data["double_well_splitting"] = split

    # Hard walls: psi' != 0 at the wall (Fourier is periodic and is left out)
    walls: dict[str, Any] = {}
    for pname, wall_sizes in HARD_WALL_SIZES.items():
        walls[pname] = {
            m: study_to_dict(stationary_study(m, pname, sizes=s)) for m, s in wall_sizes.items()
        }
    data["hard_walls"] = walls

    # Morse: a left wall too close to the classical turning region truncates the states
    morse = Morse()
    box_rows = []
    for left in (morse.box[0], morse.box[0] + 0.5):
        errs = []
        for n in (96, 192):
            e = SOLVERS["chebyshev"].solve(morse, left, morse.box[1], n, 4).energies
            errs.append(float(np.max(np.abs(e - morse.exact_energies(4)))))
        box_rows.append({"left": left, "right": morse.box[1], "errors": errs})
    data["morse_box"] = {"sizes": [96, 192], "rows": box_rows}

    # Where the fd4 time goes: LAPACK's banded driver with and without eigenvectors
    cost = []
    for n in (320, 640, 1280):
        x, h, ab = fd4_banded_matrix(ho, *ho.box, n)
        d2, e2 = 1.0 / h**2 + ho(x), np.full(n - 1, -0.5 / h**2)  # the fd2 matrix

        select = {"select": "i", "select_range": (0, 3)}
        cost.append(
            {
                "n": n,
                "fd2_eigenpairs": _clock(partial(eigh_tridiagonal, d2, e2, **select), repeats),
                "fd4_eigenvalues": _clock(
                    partial(eig_banded, ab, lower=True, eigvals_only=True, **select), repeats
                ),
                "fd4_eigenpairs": _clock(partial(eig_banded, ab, lower=True, **select), repeats),
            }
        )
    data["fd4_cost"] = cost
    return data


# --------------------------------------------------------------------------------------------
# 2. Time-dependent problem
# --------------------------------------------------------------------------------------------


def timedep_section(fig_dir: Path, quick: bool, repeats: int) -> dict[str, Any]:
    data: dict[str, Any] = {"temporal_1d": {}, "temporal_2d": {}, "spatial_1d": {}}
    fig, axes = pl.new_figure(ncols=2, width=5.6, height=4.2)
    for m in TIME_SOLVERS:
        st1 = temporal_study(m, dim=1)
        st2 = temporal_study(m, dim=2, steps=TIME_STEPS[m] if m != "rk4" else [64, 128, 256, 512])
        sp_ = spatial_study(m)
        data["temporal_1d"][m] = study_to_dict(st1)
        data["temporal_2d"][m] = study_to_dict(st2)
        data["spatial_1d"][m] = study_to_dict(sp_)
        axes[0].loglog(st1.steps, st1.errors, label=f"{m} 1D", **pl.style(m))
        axes[0].loglog(st2.steps, st2.errors, ls="--", alpha=0.7, label=f"{m} 2D", **pl.style(m))
        axes[1].loglog(sp_.steps, np.maximum(sp_.errors, 1e-16), label=m, **pl.style(m))
    pl.loglog_slope_guide(axes[0], np.array([6e-3, 0.2]), 1.5e-3, 2, "dt^2", anchor="max")
    pl.loglog_slope_guide(axes[0], np.array([2e-3, 3e-2]), 3e-7, 4, "dt^4", anchor="max")
    axes[0].set_title("Temporal error (vs exact propagator of same H)")
    axes[0].set_xlabel("time step dt")
    axes[0].set_ylabel("L2 error at t = 2")
    axes[0].legend(ncol=2)
    pl.loglog_slope_guide(axes[1], np.array([0.03, 0.5]), 3e-3, 2, "h^2", anchor="max")
    axes[1].set_title("Spatial error (vs analytic coherent state)")
    axes[1].set_xlabel("grid spacing h")
    axes[1].set_ylabel("L2 error at t = 1")
    axes[1].legend()
    pl.save(fig, fig_dir / "timedep_convergence.png")

    # Conservation (same dt for all)
    n_steps = 1000 if quick else 5000
    data["conservation_1d"] = conservation_comparison(1, 10.0, n_steps, repeats=repeats)
    data["conservation_2d"] = conservation_comparison(
        2, 5.0, 500 if quick else 1000, repeats=repeats
    )

    rho, stab = rk4_stability_demo(t_final=20.0)
    data["rk4_stability"] = {"spectral_radius": rho, "t_final": 20.0, "runs": stab}

    # Norm / energy versus time and the cost-vs-error plot (1D, N = 512, L = 30)
    grid = Grid.line(512, 30.0)
    (x,) = grid.axes()
    v = 0.5 * x**2
    psi0 = coherent_state_nd(grid, 0.0, (1.5,), (0.5,), (1.0,))
    exact = coherent_state_nd(grid, 10.0, (1.5,), (0.5,), (1.0,))
    fig, axes = pl.new_figure(ncols=3, width=5.0, height=4.0)
    cost: dict[str, Any] = {}
    for m, cls in TIME_SOLVERS.items():
        res = cls(grid, v).propagate(psi0, 10.0, n_steps)
        axes[0].semilogy(
            res.times,
            np.maximum(np.abs(res.norms - res.norms[0]), 1e-17),
            label=m,
            color=pl.STYLE[m][0],
            lw=1.5,
        )
        axes[1].semilogy(
            res.times,
            np.maximum(np.abs(res.energies - res.energies[0]) / abs(res.energies[0]), 1e-17),
            label=m,
            color=pl.STYLE[m][0],
            lw=1.5,
        )
        ks = [100, 200, 500, 1000, 2500, 5000, 10000] if m != "rk4" else [2500, 5000, 10000, 20000]
        if quick:
            ks = ks[:3]
        errs, times = [], []
        for k in ks:
            solver = cls(grid, v)
            runs = repeat_runs(partial(solver.propagate, psi0, 10.0, k, n_records=2), repeats)
            errs.append(grid.norm(runs[-1].psi - exact))
            times.append(time_stats([r.wall_time + r.setup_time for r in runs]))
        cost[m] = {"n_steps": ks, "error": errs, "seconds": times}
        axes[2].loglog([s["median"] for s in times], errs, label=m, **pl.style(m))
    data["cost_vs_error_1d"] = cost
    axes[0].set_title(f"Norm drift, N=512, dt={10.0 / n_steps:.0e}")
    axes[0].set_xlabel("t")
    axes[0].set_ylabel("abs(norm(t) - norm(0))")
    axes[1].set_title("Relative energy drift")
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("abs(E - E0) / E0")
    axes[2].set_title("Cost versus accuracy (t = 10, N = 512, L = 30)")
    axes[2].set_xlabel(f"wall time (s), median of {repeats} runs")
    axes[2].set_ylabel("L2 error vs analytic")
    for ax in axes:
        ax.legend()
    pl.save(fig, fig_dir / "conservation.png")

    data["double_well_dynamics"] = double_well_dynamics(fig_dir, quick)
    data["scattering_2d"] = scattering_2d(fig_dir, quick, repeats)
    return data


def _left_probability(
    grid: Grid, v: np.ndarray, psi0: np.ndarray, dt: float, t_end: float, n_frames: int = 0
) -> tuple[np.ndarray, np.ndarray, float, list[np.ndarray]]:
    """Split-step run recording ``P(x < 0)`` at every step (and ``n_frames`` density frames)."""
    (x,) = grid.axes()
    left = x < 0
    solver = SplitStepFourier(grid, v)
    solver.prepare(dt)
    n_steps = int(t_end / dt)
    every = max(1, n_steps // n_frames) if n_frames else 0
    psi = psi0.copy()
    norm0 = grid.norm(psi)
    p_left = np.empty(n_steps + 1)
    frames: list[np.ndarray] = []
    drift = 0.0
    for k in range(n_steps + 1):
        prob = np.abs(psi) ** 2
        p_left[k] = np.sum(prob[left]) * grid.h[0]
        if every and k % every == 0:
            frames.append(prob)
            drift = max(drift, abs(grid.norm(psi) - norm0))
        psi = solver.step(psi)
    drift = max(drift, abs(grid.norm(psi) - norm0))
    return dt * np.arange(n_steps + 1), p_left, drift, frames


def _fit_frequency(
    t: np.ndarray, y: np.ndarray, guess: float, filter_width: float
) -> tuple[float, np.ndarray, np.ndarray]:
    """Frequency of the slow oscillation in ``y(t)``.

    A Gaussian low-pass filter of width ``filter_width`` removes the fast intra-well
    oscillations; ``a + b cos(w t) + c sin(w t)`` is then fitted by least squares (linear in
    ``a, b, c``, a 1D minimisation in ``w``) away from the filtered edges.
    """
    dt = float(t[1] - t[0])
    smooth = gaussian_filter1d(y, filter_width / dt, mode="nearest")
    cut = int(4 * filter_width / dt)
    tt, yy = t[cut:-cut], smooth[cut:-cut]

    def residual(w: float) -> float:
        a = np.column_stack([np.ones_like(tt), np.cos(w * tt), np.sin(w * tt)])
        coef = np.linalg.lstsq(a, yy, rcond=None)[0]
        return float(np.sum((a @ coef - yy) ** 2))

    w = minimize_scalar(residual, bracket=(0.95 * guess, guess, 1.05 * guess), tol=1e-14).x
    return float(w), tt, yy


def _propagator_splitting(grid: Grid, v: np.ndarray, dt: float, states: np.ndarray) -> float:
    """``E1 - E0`` of the one-step Strang propagator from its eigenphases (dense, small grids).

    The two eigenvectors with the largest overlap with ``states[:, 0]`` and ``states[:, 1]``
    are used, so the ``2 pi / dt`` ambiguity of the phases cannot pick a wrong pair.
    """
    solver = SplitStepFourier(grid, v)
    solver.prepare(dt)
    eye = np.eye(grid.size, dtype=complex)
    u = np.column_stack([solver.step(eye[:, j]) for j in range(grid.size)])
    lam, vecs = np.linalg.eig(u)
    energies = []
    for k in range(2):
        j = int(np.argmax(np.abs(vecs.conj().T @ states[:, k])))
        energies.append(-np.angle(lam[j]) / dt)
    return float(energies[1] - energies[0])


def double_well_dynamics(fig_dir: Path, quick: bool) -> dict[str, Any]:
    """Tunnel splitting from the dynamics of a Gaussian packet, against the Chebyshev solver.

    The initial state is a Gaussian at the left minimum with the width of the local harmonic
    ground state, built without any eigensolver.  The slow oscillation of ``P(x < 0)`` is
    fitted (see :func:`_fit_frequency`), and its frequency is compared with ``E1 - E0`` from the
    Chebyshev solver, a different basis and method.  For each ``dt`` the splitting of the Strang
    propagator itself, from the eigenphases of its one-step matrix, is listed as well: the fit
    should reproduce it, and both should differ from the stationary value by ``O(dt^2)``.
    """
    dw = DoubleWell()
    n, length = 256, 12.0
    grid = Grid.line(n, length)
    (x,) = grid.axes()
    v = dw(x)
    ref = SOLVERS["chebyshev"].solve(dw, -length / 2, length / 2, 200, 2).energies
    split_ref = float(ref[1] - ref[0])
    states = SOLVERS["fourier"].solve(dw, -length / 2, length / 2, n, 2).states.astype(complex)
    omega = np.sqrt(8.0 * dw.lam * dw.a**2)  # V''(-a)
    psi0 = grid.normalize(np.exp(-0.5 * omega * (x + dw.a) ** 2).astype(complex))
    filter_width, periods = 15.0, 1.1
    t_end = periods * 2.0 * np.pi / split_ref
    rows: list[dict[str, float | None]] = []
    plot = None
    for dt in (0.04, 0.02) if quick else (0.04, 0.02, 0.01):
        t, p_left, drift, frames = _left_probability(
            grid, v, psi0, dt, t_end, n_frames=400 if dt == 0.02 else 0
        )
        w, tt, yy = _fit_frequency(t, p_left, split_ref, filter_width)
        prop = _propagator_splitting(grid, v, dt, states)
        rows.append(
            {
                "dt": dt,
                "splitting_fit": w,
                "relative_difference_fit": (w - split_ref) / split_ref,
                "splitting_propagator": prop,
                "relative_difference_propagator": (prop - split_ref) / split_ref,
                "max_norm_drift": drift,
            }
        )
        if frames:
            plot = (t, p_left, tt, yy, frames)
    rel = [(float(r["dt"] or 0.0), float(r["relative_difference_fit"] or 0.0)) for r in rows]
    rows[0]["order"] = None
    for i in range(1, len(rows)):
        (dt_a, a), (dt_b, b) = rel[i - 1], rel[i]
        rows[i]["order"] = float(np.log(a / b) / np.log(dt_a / dt_b))

    assert plot is not None
    t, p_left, tt, yy, frames = plot
    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    ax = axes[0]
    ax.imshow(np.array(frames), extent=(x[0], x[-1], t[-1], t[0]), aspect="auto", cmap="viridis")
    ax.set_xlim(-4.5, 4.5)
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title("Density |psi(x,t)|^2, Gaussian start in the left well")
    ax.grid(False)
    ax = axes[1]
    ax.plot(t, p_left, color=pl.STYLE["crank-nicolson"][0], lw=0.4, alpha=0.6, label="P(x<0)")
    ax.plot(tt, yy, color=pl.INK, lw=1.4, label=f"low-pass filtered ({filter_width:g} time units)")
    ax.axvline(
        np.pi / split_ref,
        color=pl.STYLE["split-step"][0],
        ls="--",
        lw=1.2,
        label="pi/(E1-E0) from the Chebyshev solver",
    )
    ax.set_xlabel("t")
    ax.set_ylabel("probability in left well")
    ax.set_title("Tunnelling oscillation (split-step, dt = 0.02)")
    ax.legend(loc="upper center")
    pl.save(fig, fig_dir / "double_well_tunnelling.png")
    return {
        "n": n,
        "length": length,
        "reference": "chebyshev, N = 200",
        "splitting_reference": split_ref,
        "filter_width": filter_width,
        "periods": periods,
        "rows": rows,
    }


def scattering_2d(fig_dir: Path, quick: bool, repeats: int) -> dict[str, Any]:
    """2D packet hitting a Gaussian bump: split-step versus Crank-Nicolson."""
    n, length = (64, 24.0) if quick else (128, 24.0)
    grid = Grid.square(n, length)
    xx, yy = grid.mesh()
    v = 2.5 * np.exp(-(xx**2 + yy**2) / 2.0)
    psi0 = gaussian_packet(grid, (-6.0, 0.0), 1.2, (2.0, 0.0))
    t_final, dt = 4.0, 0.01
    n_steps = int(t_final / dt)
    results = {}
    seconds = {}
    for cls in (SplitStepFourier, CrankNicolson):
        s = cls(grid, v)
        runs = repeat_runs(partial(s.propagate, psi0, t_final, n_steps, n_records=20), repeats)
        results[s.name] = runs[-1]
        seconds[s.name] = time_stats([r.wall_time + r.setup_time for r in runs])
    ss, cn = results["split-step"], results["crank-nicolson"]
    trans = {
        k: float(np.sum(np.abs(r.psi[xx > 0]) ** 2) * grid.cell_volume) for k, r in results.items()
    }
    rel = float(grid.norm(ss.psi - cn.psi) / grid.norm(ss.psi))
    rel_coarse = float("nan")
    if not quick:  # same comparison with the grid spacing doubled
        g2 = Grid.square(n // 2, length)
        x2, y2 = g2.mesh()
        v2 = 2.5 * np.exp(-(x2**2 + y2**2) / 2.0)
        p2 = gaussian_packet(g2, (-6.0, 0.0), 1.2, (2.0, 0.0))
        a = SplitStepFourier(g2, v2).propagate(p2, t_final, n_steps, n_records=2).psi
        b = CrankNicolson(g2, v2).propagate(p2, t_final, n_steps, n_records=2).psi
        rel_coarse = float(g2.norm(a - b) / g2.norm(a))
    fig, axes = pl.new_figure(ncols=3, width=4.2, height=4.0)
    snaps = {}
    s = SplitStepFourier(grid, v)
    s.prepare(dt)
    psi = psi0.copy()
    done = 0
    for tt in (0.0, 2.0, 4.0):
        target = round(tt / dt)
        for _ in range(target - done):
            psi = s.step(psi)
        done = target
        snaps[tt] = np.abs(psi) ** 2
    for ax, (tt, d) in zip(axes, snaps.items(), strict=True):
        ax.imshow(
            d.T,
            origin="lower",
            extent=(-length / 2, length / 2, -length / 2, length / 2),
            cmap="viridis",
        )
        ax.contour(xx.T, yy.T, v.T, levels=[1.25], colors="white", linewidths=0.8)
        ax.set_xlim(-10, 10)
        ax.set_ylim(-10, 10)
        ax.set_title(f"t = {tt:.0f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(False)
    pl.save(fig, fig_dir / "wavepacket_2d.png")
    return {
        "n": n,
        "length": length,
        "dt": dt,
        "n_steps": n_steps,
        "transmitted_split_step": trans["split-step"],
        "transmitted_crank_nicolson": trans["crank-nicolson"],
        "relative_l2_difference": rel,
        "relative_l2_difference_half_grid": rel_coarse,
        "norm_drift_split_step": ss.max_norm_drift,
        "norm_drift_crank_nicolson": cn.max_norm_drift,
        "seconds_split_step": seconds["split-step"],
        "seconds_crank_nicolson": seconds["crank-nicolson"],
    }


# --------------------------------------------------------------------------------------------
# 3. Tunnelling
# --------------------------------------------------------------------------------------------


def tunneling_section(fig_dir: Path, quick: bool, repeats: int) -> dict[str, Any]:
    barrier = RectangularBarrier(height=1.0, width=1.55)
    energies = np.linspace(0.05, 3.0, 300)
    # FD curve at m = 81 points per barrier width
    m_curve = 81
    h = barrier.width / m_curve
    curve = transmission_fd(energies, barrier.height, m_curve, h)
    exact_curve = barrier.transmission(energies)

    ms = [5, 11, 21, 41, 81, 161]
    e_probe = 0.72
    rows = []
    hs, errs = [], []
    for m in ms:
        hh = barrier.width / m
        r = transmission_fd(e_probe, barrier.height, m, hh)
        ex = float(barrier.transmission(e_probe)[0])
        hs.append(hh)
        errs.append(abs(float(r.transmission[0]) - ex))
        rows.append(
            (
                m,
                hh,
                float(r.transmission[0]),
                ex,
                errs[-1],
                float(r.reflection[0] + r.transmission[0] - 1),
            )
        )
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(
        np.array(hs[:-1]) / np.array(hs[1:])
    )
    k0s = [0.6, 0.8, 1.0, 1.2, 1.6, 2.0] if not quick else [0.8, 1.2]
    pk = [packet_transmission(k0, barrier.height, barrier.width) for k0 in k0s]
    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    ax = axes[0]
    ax.plot(energies, exact_curve, color=pl.INK, lw=1.6, label="exact T(E)")
    ax.plot(
        energies[::6],
        curve.transmission[::6],
        color=pl.STYLE["fd2"][0],
        lw=0,
        marker="o",
        markersize=3.5,
        label=f"FD scattering ({m_curve} pts in barrier)",
    )
    ax.plot(
        [0.5 * p.k0**2 for p in pk],
        [p.numeric for p in pk],
        color=pl.STYLE["split-step"][0],
        lw=0,
        marker="s",
        markersize=6,
        label="wave packet (split-step)",
    )
    ax.axvline(barrier.height, color="#898781", ls=":", lw=1)
    ax.set_xlabel("E (V0 = 1, a = 1.55)")
    ax.set_ylabel("transmission T")
    ax.set_title("Rectangular barrier")
    ax.legend(loc="upper center")
    ax = axes[1]
    ax.loglog(hs, errs, label="FD scattering, E = 0.72", **pl.style("fd2"))
    h_range = np.array([min(hs), max(hs)])
    pl.loglog_slope_guide(ax, h_range, errs[-1] * 0.4, 2, "h^2")
    ax.set_xlabel("grid spacing h")
    ax.set_ylabel("abs(T_FD - T_exact)")
    ax.set_title("Convergence of the numerical T")
    ax.legend()
    pl.save(fig, fig_dir / "tunneling.png")
    return {
        "energy": e_probe,
        "fd_convergence": [
            dict(
                m=r[0],
                h=r[1],
                T_fd=r[2],
                T_exact=r[3],
                error=r[4],
                order=float(orders[i - 1]) if i else None,
                flux_defect=r[5],
            )
            for i, r in enumerate(rows)
        ],
        "fd_fitted_order": float(np.polyfit(np.log(hs), np.log(errs), 1)[0]),
        "packet": [p.__dict__ for p in pk],
        "packet_later": [packet_transmission_later(k0) for k0 in ([0.6] if quick else [0.6, 1.2])],
    }


# --------------------------------------------------------------------------------------------
# 4. Poisson / heat
# --------------------------------------------------------------------------------------------


def elliptic_section(fig_dir: Path, quick: bool, repeats: int) -> dict[str, Any]:
    sizes = [31, 63, 127, 255] if not quick else [31, 63]
    methods = POISSON_SOLVERS
    iters: dict[str, list[int]] = {m: [] for m in methods}
    secs: dict[str, list[dict[str, float]]] = {m: [] for m in methods}
    disc_err = []
    hist: dict[tuple[str, int], list[float]] = {}
    for n in sizes:
        f, u = manufactured_problem(n)
        for m in methods:
            runs = repeat_runs(partial(solve_poisson, f, m, tol=1e-10), repeats)
            r = runs[-1]
            iters[m].append(r.iterations)
            secs[m].append(time_stats([q.setup_time + q.solve_time for q in runs]))
            if n == sizes[-1] or n == 127:
                hist[(m, n)] = r.residuals
            if m == "direct":
                disc_err.append(l2_error(r.x, u))
    hs = np.array([1.0 / (n + 1) for n in sizes])
    d_orders = np.log(np.array(disc_err[:-1]) / np.array(disc_err[1:])) / np.log(hs[:-1] / hs[1:])
    data: dict[str, Any] = {
        "poisson": {
            "sizes": sizes,
            "tolerance": 1e-10,
            "iterations": iters,
            "seconds": secs,
            "discretisation_error": disc_err,
            "discretisation_orders": d_orders.tolist(),
        }
    }

    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    n_hist = 127 if 127 in sizes else sizes[-1]
    for m in ("cg", "cg-ic0", "mg", "mg-cg"):
        res = hist[(m, n_hist)]
        axes[0].semilogy(np.arange(len(res)), res, label=m, color=pl.STYLE[m][0], lw=1.6)
    axes[0].set_xlabel("iteration")
    axes[0].set_ylabel("relative residual")
    axes[0].set_title(f"Residual history, n = {n_hist}")
    axes[0].set_xlim(0, 120)
    axes[0].legend()
    for m in ("cg", "cg-ic0", "mg", "mg-cg"):
        axes[1].loglog(sizes, iters[m], label=m, **pl.style(m))
    axes[1].set_xlabel("n (grid is n x n)")
    axes[1].set_ylabel("iterations to relative residual 1e-10")
    axes[1].set_title("Iteration count versus grid size")
    axes[1].legend()
    pl.save(fig, fig_dir / "poisson_solvers.png")

    # Heat equation: temporal order
    n = 31
    t_final = 0.1
    ks = [10, 20, 40, 80]
    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    heat_orders = {}
    for scheme in ("backward-euler", "crank-nicolson"):
        u0, ex = heat_exact_semidiscrete(n, (1, 2), t_final)
        errs = []
        for k in ks:
            hr = heat_solve(u0, t_final, k, scheme, "direct")
            errs.append(float(np.linalg.norm(hr.u - ex) / (n + 1)))
        dts = [t_final / k for k in ks]
        axes[0].loglog(dts, errs, label=scheme, **pl.style(scheme))
        o = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(
            np.array(dts[:-1]) / np.array(dts[1:])
        )
        heat_orders[scheme] = {"dt": dts, "error": errs, "orders": o.tolist()}
    axes[0].set_xlabel("time step dt")
    axes[0].set_ylabel("L2 error at t = 0.1")
    axes[0].set_title("Heat equation: temporal convergence")
    dd = np.array([t_final / 80, t_final / 10])
    pl.loglog_slope_guide(axes[0], dd, 2e-3, 1, "dt", anchor="max")
    pl.loglog_slope_guide(axes[0], dd, 1e-4, 2, "dt^2", anchor="max")
    axes[0].legend()

    # Linear solvers inside the CN time loop, each set up once per run
    n2 = 63 if not quick else 31
    u_init = np.zeros((n2, n2))
    u_init[n2 // 3 : n2 // 3 + 6, n2 // 2 : n2 // 2 + 9] = 1.0
    n_heat_steps = 20
    rows = []
    sols = {}
    for ls in ("direct", "cg", "cg-ic0", "mg-cg"):
        heat_runs = repeat_runs(
            partial(heat_solve, u_init, 0.02, n_heat_steps, "crank-nicolson", ls, tol=1e-12),
            repeats,
        )
        hr = heat_runs[-1]
        sols[ls] = hr.u
        rows.append(
            {
                "method": ls,
                "iterations": hr.total_inner_iterations,
                "setup_seconds": time_stats([q.setup_time for q in heat_runs]),
                "seconds": time_stats([q.wall_time for q in heat_runs]),
                "max_difference_to_direct": float(np.max(np.abs(hr.u - sols["direct"]))),
            }
        )
    u_final = sols["mg-cg"]
    ax = axes[1]
    ax.imshow(u_final.T, origin="lower", extent=(0, 1, 0, 1), cmap="viridis")
    ax.contour(
        u_init.T, levels=[0.5], origin="lower", extent=(0, 1, 0, 1), colors="white", linewidths=0.8
    )
    ax.set_title("Hot patch at t = 0.02 (white: initial edge)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    pl.save(fig, fig_dir / "heat_equation.png")
    data["heat"] = {
        "n": n,
        "t_final": t_final,
        "n_steps": ks,
        "orders": heat_orders,
        "solvers": {"n": n2, "n_steps": n_heat_steps, "rows": rows},
    }
    return data


# --------------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------------


def generate_report(
    out_dir: Path,
    readme: Path | None = None,
    quick: bool = False,
    verbose: bool = True,
    repeats: int = REPEATS,
) -> dict[str, Any]:
    """Run all studies; write ``figures/``, ``results.json``, ``tables.md`` under ``out_dir``."""
    out_dir = Path(out_dir)
    fig_dir = out_dir / "figures"
    res_dir = out_dir / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {
        "environment": environment(),
        "quick": quick,
        "timing": {"repeats": repeats, "statistic": "median, min and max after one warm-up run"},
    }
    sections = {
        "stationary": stationary_section,
        "timedep": timedep_section,
        "tunneling": tunneling_section,
        "elliptic": elliptic_section,
    }
    for name, fn in sections.items():
        t0 = time.perf_counter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            results[name] = fn(fig_dir, quick, repeats)
        if verbose:
            print(f"[report] {name}: {time.perf_counter() - t0:.1f} s")
    results = _clean(results)
    (res_dir / "results.json").write_text(
        json.dumps(results, indent=1, allow_nan=False) + "\n", encoding="utf-8"
    )
    tables = render_tables(results)
    (res_dir / "tables.md").write_text(tables_markdown(tables), encoding="utf-8")
    if readme is not None:
        found = update_readme(Path(readme), tables)
        if verbose:
            print(f"[report] README sections updated: {', '.join(found)}")
            missing = sorted(set(tables) - set(found))
            if missing:
                print(f"[report] no README marker for: {', '.join(missing)}")
    return results
