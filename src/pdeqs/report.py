"""Run every study, write PNG figures, JSON results and markdown tables (``pdeqs report``).

Nothing here is hand-typed: tables and figures are produced from the runs, and the README tables
are refreshed between ``<!-- BEGIN:name -->`` / ``<!-- END:name -->`` markers.
"""

from __future__ import annotations

import json
import os
import platform
import re
import time
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import scipy

from . import plotting as pl
from .convergence import (
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
from .potentials import DoubleWell, Harmonic, RectangularBarrier
from .stationary import SOLVERS, align_sign
from .timedep import (
    TIME_SOLVERS,
    CrankNicolson,
    SplitStepFourier,
    coherent_state_nd,
    gaussian_packet,
)
from .tunneling import packet_transmission, transmission_fd

# --------------------------------------------------------------------------------------------
# Formatting helpers
# --------------------------------------------------------------------------------------------


def sci(x: float, digits: int = 2) -> str:
    return "-" if x != x else f"{x:.{digits}e}"


def num(x: float, digits: int = 2) -> str:
    return "-" if x != x else f"{x:.{digits}f}"


def md_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def environment() -> dict[str, str]:
    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "cpu_count": str(os.cpu_count()),
    }


def update_readme(readme: Path, tables: dict[str, str]) -> list[str]:
    """Replace the text between BEGIN/END markers; returns the names that were found."""
    text = readme.read_text(encoding="utf-8")
    found = []
    for name, body in tables.items():
        pattern = re.compile(
            rf"(<!-- BEGIN:{re.escape(name)} -->)(.*?)(<!-- END:{re.escape(name)} -->)", re.S
        )
        if pattern.search(text):
            text = pattern.sub(lambda m, b=body: f"{m.group(1)}\n{b}\n{m.group(3)}", text)  # type: ignore[misc]
            found.append(name)
    readme.write_text(text, encoding="utf-8")
    return found


def study_to_dict(s: ConvergenceStudy) -> dict[str, Any]:
    return {
        "name": s.name,
        "parameter": s.parameter,
        "reference": s.reference,
        "expected_order": s.expected,
        "fitted_order": s.fitted_order,
        "rows": s.rows(),
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


def stationary_section(fig_dir: Path, quick: bool) -> tuple[dict[str, Any], dict[str, str]]:
    data: dict[str, Any] = {}
    tables: dict[str, str] = {}
    fig, axes = pl.new_figure(ncols=3, width=5.0, height=4.0)
    for ax, pname in zip(axes, POTENTIALS, strict=True):
        data[pname] = {}
        rows = []
        for method in SOLVERS:
            sizes = STATIONARY_SIZES[method][:4] if quick else None
            st = stationary_study(method, pname, sizes=sizes)
            data[pname][method] = study_to_dict(st)
            n_arr = np.array(st.sizes)
            err = np.maximum(st.errors, 1e-16)
            ax.loglog(n_arr, err, label=method, **pl.style(method))
            last = st.pairwise_orders[-1] if len(st.errors) > 1 else float("nan")
            rows.append(
                [
                    method,
                    "spectral" if st.is_spectral else st.expected,
                    "-" if st.is_spectral else num(st.fitted_order),
                    "-" if st.is_spectral else num(last),
                    f"{sci(st.errors[0])} (N={st.sizes[0]})",
                    f"{sci(st.errors[-1])} (N={st.sizes[-1]})",
                    num(st.wall_times[-1] * 1e3, 1),
                ]
            )
        ref_name = data[pname]["fd2"]["reference"]
        ax.set_title(f"{pname}  (reference: {ref_name})")
        ax.set_xlabel("N (interior nodes / modes)")
        ax.set_ylabel("max error of lowest 4 levels")
        ax.set_ylim(1e-15, 1e1)
        if pname == "harmonic":
            n_ref = np.array([40.0, 1280.0])
            pl.loglog_slope_guide(ax, n_ref, 3e-1, -2, "N^-2")
            pl.loglog_slope_guide(ax, n_ref, 5e-3, -4, "N^-4")
        ax.legend(loc="lower left", ncol=2)
        tables[f"stationary_{pname.replace('-', '_')}"] = md_table(
            [
                "Method",
                "Nominal order",
                "Fitted order",
                "Last-pair order",
                "Error, smallest N",
                "Error, largest N",
                "Time, largest N (ms)",
            ],
            rows,
        )
    pl.save(fig, fig_dir / "stationary_convergence.png")

    # Wavefunction accuracy against the exact oscillator eigenfunctions, plus eigenfunction plot
    ho = Harmonic()
    rows = []
    wf: dict[str, Any] = {}
    for method, solver in SOLVERS.items():
        res = solver.solve(ho, *ho.box, 200, 6)
        errs = []
        for k in (0, 3):
            ex = ho.exact_state(res.x, k)
            psi = align_sign(res.states[:, k], ex)
            errs.append(float(np.max(np.abs(psi - ex))))
        wf[method] = errs
        rows.append(
            [
                method,
                sci(errs[0]),
                sci(errs[1]),
                sci(float(np.max(np.abs(res.energies - ho.exact_energies(6))))),
            ]
        )
    data["harmonic_wavefunctions_N200"] = wf
    tables["stationary_wavefunctions"] = md_table(
        ["Method", "max error of psi_0", "max error of psi_3", "max energy error (6 levels)"],
        rows,
    )

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
    rows = []
    split: dict[str, float] = {}
    for method in ("fd2", "fd4", "numerov", "shooting", "fourier", "chebyshev"):
        n = 400 if method in ("fd2", "fd4", "numerov", "shooting") else 96
        r = SOLVERS[method].solve(dw, *dw.box, n, 2)
        split[method] = float(r.energies[1] - r.energies[0])
        rows.append([method, str(n), f"{r.energies[0]:.8f}", f"{split[method]:.6e}"])
    data["double_well_splitting"] = split
    tables["double_well_splitting"] = md_table(
        ["Method", "N", "E_0", "E_1 - E_0 (tunnel splitting)"], rows
    )
    return data, tables


# --------------------------------------------------------------------------------------------
# 2. Time-dependent problem
# --------------------------------------------------------------------------------------------


def timedep_section(fig_dir: Path, quick: bool) -> tuple[dict[str, Any], dict[str, str]]:
    data: dict[str, Any] = {"temporal_1d": {}, "temporal_2d": {}, "spatial_1d": {}}
    tables: dict[str, str] = {}
    methods = list(TIME_SOLVERS)
    fig, axes = pl.new_figure(ncols=2, width=5.6, height=4.2)
    order_rows = []
    for m in methods:
        st1 = temporal_study(m, dim=1)
        st2 = temporal_study(m, dim=2, steps=TIME_STEPS[m] if m != "rk4" else [64, 128, 256, 512])
        sp_ = spatial_study(m)
        data["temporal_1d"][m] = study_to_dict(st1)
        data["temporal_2d"][m] = study_to_dict(st2)
        data["spatial_1d"][m] = study_to_dict(sp_)
        axes[0].loglog(st1.steps, st1.errors, label=f"{m} 1D", **pl.style(m))
        axes[0].loglog(st2.steps, st2.errors, ls="--", alpha=0.7, label=f"{m} 2D", **pl.style(m))
        axes[1].loglog(sp_.steps, np.maximum(sp_.errors, 1e-16), label=m, **pl.style(m))
        order_rows.append(
            [
                m,
                str(TIME_SOLVERS[m].order),
                num(st1.fitted_order),
                num(st1.pairwise_orders[-1]),
                num(st2.fitted_order),
                num(st2.pairwise_orders[-1]),
                "spectral" if sp_.is_spectral else num(sp_.fitted_order),
                "-" if sp_.is_spectral else num(sp_.pairwise_orders[-1]),
            ]
        )
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
    tables["timedep_orders"] = md_table(
        [
            "Method",
            "Nominal dt-order",
            "1D fitted",
            "1D last pair",
            "2D fitted",
            "2D last pair",
            "h-order fitted (1D)",
            "h-order last pair",
        ],
        order_rows,
    )

    # Conservation (same dt for all)
    n_steps = 1000 if quick else 5000
    cons1 = conservation_comparison(1, 10.0, n_steps)
    cons2 = conservation_comparison(2, 5.0, 500 if quick else 1000)
    data["conservation_1d"] = cons1
    data["conservation_2d"] = cons2

    def cons_table(rows: list[dict[str, Any]]) -> str:
        return md_table(
            [
                "Method",
                "dt",
                "max abs(norm - 1)",
                "max rel. energy drift",
                "L2 error vs analytic",
                "Time (s)",
            ],
            [
                [
                    r["method"],
                    sci(r["dt"], 1),
                    sci(r["max_norm_drift"]),
                    sci(r["relative_energy_drift"]),
                    sci(r["l2_error_vs_analytic"]),
                    num(r["seconds"], 3),
                ]
                for r in rows
            ],
        )

    tables["conservation_1d"] = cons_table(cons1)
    tables["conservation_2d"] = cons_table(cons2)

    rho, stab = rk4_stability_demo()
    data["rk4_stability"] = {"spectral_radius": rho, "runs": stab}
    tables["rk4_stability"] = md_table(
        ["dt", "dt x rho(H)", "norm after t = 2"],
        [[sci(r["dt"], 2), num(r["dt_times_rho"], 2), sci(r["final_norm"], 3)] for r in stab],
    )

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
            np.maximum(np.abs(res.norms - 1.0), 1e-17),
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
            r = cls(grid, v).propagate(psi0, 10.0, k, n_records=2)
            errs.append(grid.norm(r.psi - exact))
            times.append(r.wall_time + r.setup_time)
        cost[m] = {"n_steps": ks, "error": errs, "seconds": times}
        axes[2].loglog(times, errs, label=m, **pl.style(m))
    data["cost_vs_error_1d"] = cost
    axes[0].set_title(f"Norm drift, N=512, dt={10.0 / n_steps:.0e}")
    axes[0].set_xlabel("t")
    axes[0].set_ylabel("abs(norm - 1)")
    axes[1].set_title("Relative energy drift")
    axes[1].set_xlabel("t")
    axes[1].set_ylabel("abs(E - E0) / E0")
    axes[2].set_title("Cost versus accuracy (t = 10, N = 512, L = 30)")
    axes[2].set_xlabel("wall time (s)")
    axes[2].set_ylabel("L2 error vs analytic")
    for ax in axes:
        ax.legend()
    pl.save(fig, fig_dir / "conservation.png")

    # Double-well tunnelling oscillation: time-dependent vs stationary splitting
    extra = double_well_dynamics(fig_dir, quick)
    data["double_well_dynamics"] = extra
    tables["double_well_dynamics"] = md_table(
        ["Quantity", "Value"],
        [
            [
                "splitting E1-E0 from stationary Fourier solver",
                f"{extra['splitting_stationary']:.6e}",
            ],
            [
                "splitting from time-dependent first minimum (pi / t_min)",
                f"{extra['splitting_dynamic']:.6e}",
            ],
            ["relative difference", sci(extra["relative_difference"])],
            ["max abs(norm - 1) over the run", sci(extra["max_norm_drift"])],
        ],
    )

    scat = scattering_2d(fig_dir, quick)
    data["scattering_2d"] = scat
    tables["scattering_2d"] = md_table(
        ["Quantity", "Value"],
        [
            ["grid", f"{scat['n']} x {scat['n']}, L = {scat['length']}"],
            ["dt, steps", f"{scat['dt']}, {scat['n_steps']}"],
            ["transmitted probability (split-step)", f"{scat['transmitted_split_step']:.5f}"],
            [
                "transmitted probability (Crank-Nicolson)",
                f"{scat['transmitted_crank_nicolson']:.5f}",
            ],
            [
                "relative L2 difference of the two wave functions",
                sci(scat["relative_l2_difference"]),
            ],
            [
                "same difference with the grid spacing doubled (half as many points)",
                sci(scat["relative_l2_difference_half_grid"]),
            ],
            [
                "max abs(norm - 1): split-step / CN",
                f"{sci(scat['norm_drift_split_step'])} / {sci(scat['norm_drift_crank_nicolson'])}",
            ],
            [
                "wall time: split-step / CN (s)",
                f"{scat['seconds_split_step']:.2f} / {scat['seconds_crank_nicolson']:.2f}",
            ],
        ],
    )
    return data, tables


def double_well_dynamics(fig_dir: Path, quick: bool) -> dict[str, float]:
    """Packet in the left well of the double well; compares the oscillation with E1 - E0."""
    dw = DoubleWell()
    n, length = 256, 12.0
    grid = Grid.line(n, length)
    (x,) = grid.axes()
    ref = SOLVERS["fourier"].solve(dw, -length / 2, length / 2, n, 2)
    e0, e1 = ref.energies
    p0, p1 = ref.states[:, 0], ref.states[:, 1]
    # left-localised combination; sign chosen so that most weight sits at x < 0
    psi0 = (p0 + p1) / np.sqrt(2.0)
    if np.sum(psi0[x < 0] ** 2) < np.sum(psi0[x > 0] ** 2):
        psi0 = (p0 - p1) / np.sqrt(2.0)
    psi0 = grid.normalize(psi0.astype(complex))
    solver = SplitStepFourier(grid, dw(x))
    dt = 0.02
    t_end = 1.15 * np.pi / (e1 - e0)
    n_steps = int(t_end / dt)
    solver.prepare(dt)
    n_frames = 400
    every = max(1, n_steps // n_frames)
    psi = psi0.copy()
    times, left, dens, norms = [], [], [], []
    for k in range(n_steps + 1):
        if k % every == 0:
            times.append(k * dt)
            prob = np.abs(psi) ** 2
            left.append(float(np.sum(prob[x < 0]) * grid.h[0]))
            dens.append(prob)
            norms.append(grid.norm(psi))
        psi = solver.step(psi)
    times_a, left_a = np.array(times), np.array(left)
    i = int(np.argmin(left_a))
    # parabolic refinement of the minimum
    if 0 < i < len(left_a) - 1:
        y0, y1, y2 = left_a[i - 1], left_a[i], left_a[i + 1]
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        t_min = times_a[i] + shift * (times_a[1] - times_a[0])
    else:
        t_min = times_a[i]
    split_dyn = float(np.pi / t_min)
    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    ax = axes[0]
    ax.imshow(
        np.array(dens),
        extent=(x[0], x[-1], times_a[-1], times_a[0]),
        aspect="auto",
        cmap="viridis",
    )
    ax.set_xlim(-4.5, 4.5)
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title("Density |psi(x,t)|^2 in the double well")
    ax.grid(False)
    ax = axes[1]
    ax.plot(times_a, left_a, color=pl.STYLE["crank-nicolson"][0], lw=1.6, label="split-step P(x<0)")
    ax.axvline(
        np.pi / (e1 - e0),
        color=pl.STYLE["split-step"][0],
        ls="--",
        lw=1.2,
        label="pi/(E1-E0) from stationary solver",
    )
    ax.set_xlabel("t")
    ax.set_ylabel("probability in left well")
    ax.set_title("Tunnelling oscillation")
    ax.legend(loc="upper right")
    pl.save(fig, fig_dir / "double_well_tunnelling.png")
    return {
        "splitting_stationary": float(e1 - e0),
        "splitting_dynamic": split_dyn,
        "relative_difference": float(abs(split_dyn - (e1 - e0)) / (e1 - e0)),
        "max_norm_drift": float(np.max(np.abs(np.array(norms) - 1.0))),
        "dt": dt,
        "n_steps": float(n_steps),
    }


def scattering_2d(fig_dir: Path, quick: bool) -> dict[str, Any]:
    """2D packet hitting a Gaussian bump: split-step versus Crank-Nicolson."""
    n, length = (64, 24.0) if quick else (128, 24.0)
    grid = Grid.square(n, length)
    xx, yy = grid.mesh()
    v = 2.5 * np.exp(-(xx**2 + yy**2) / 2.0)
    psi0 = gaussian_packet(grid, (-6.0, 0.0), 1.2, (2.0, 0.0))
    t_final, dt = 4.0, 0.01
    n_steps = int(t_final / dt)
    results = {}
    for cls in (SplitStepFourier, CrankNicolson):
        s = cls(grid, v)
        results[s.name] = s.propagate(psi0, t_final, n_steps, n_records=20)
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
        "seconds_split_step": ss.wall_time + ss.setup_time,
        "seconds_crank_nicolson": cn.wall_time + cn.setup_time,
    }


# --------------------------------------------------------------------------------------------
# 3. Tunnelling
# --------------------------------------------------------------------------------------------


def tunneling_section(fig_dir: Path, quick: bool) -> tuple[dict[str, Any], dict[str, str]]:
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
    fd_table = md_table(
        ["Points in barrier", "h", "T (FD)", "T (exact)", "abs error", "order", "abs(R + T - 1)"],
        [
            [
                str(r[0]),
                sci(r[1], 2),
                f"{r[2]:.6f}",
                f"{r[3]:.6f}",
                sci(r[4]),
                "-" if i == 0 else num(float(orders[i - 1])),
                sci(abs(r[5]), 1),
            ]
            for i, r in enumerate(rows)
        ],
    )

    k0s = [0.6, 0.8, 1.0, 1.2, 1.6, 2.0] if not quick else [0.8, 1.2]
    pk = [packet_transmission(k0, barrier.height, barrier.width) for k0 in k0s]
    packet_table = md_table(
        [
            "k0",
            "mean energy k0^2/2",
            "T(mean energy)",
            "T averaged over packet (exact)",
            "T from split-step",
            "abs difference",
        ],
        [
            [
                num(p.k0, 2),
                num(0.5 * p.k0**2, 3),
                f"{p.analytic_mean_energy:.5f}",
                f"{p.analytic_averaged:.5f}",
                f"{p.numeric:.5f}",
                sci(abs(p.numeric - p.analytic_averaged), 1),
            ]
            for p in pk
        ],
    )
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
    ax.legend(loc="lower right")
    ax = axes[1]
    ax.loglog(hs, errs, label="FD scattering, E = 0.72", **pl.style("fd2"))
    h_range = np.array([min(hs), max(hs)])
    pl.loglog_slope_guide(ax, h_range, errs[-1] * 0.4, 2, "h^2")
    ax.set_xlabel("grid spacing h")
    ax.set_ylabel("abs(T_FD - T_exact)")
    ax.set_title("Convergence of the numerical T")
    ax.legend()
    pl.save(fig, fig_dir / "tunneling.png")
    data = {
        "fd_convergence": [
            dict(m=r[0], h=r[1], T_fd=r[2], T_exact=r[3], error=r[4], flux_defect=r[5])
            for r in rows
        ],
        "fd_fitted_order": float(np.polyfit(np.log(hs), np.log(errs), 1)[0]),
        "packet": [p.__dict__ for p in pk],
    }
    return data, {"tunneling_fd": fd_table, "tunneling_packet": packet_table}


# --------------------------------------------------------------------------------------------
# 4. Poisson / heat
# --------------------------------------------------------------------------------------------


def elliptic_section(fig_dir: Path, quick: bool) -> tuple[dict[str, Any], dict[str, str]]:
    sizes = [31, 63, 127, 255] if not quick else [31, 63]
    methods = POISSON_SOLVERS
    data: dict[str, Any] = {"poisson": {}, "heat": {}}
    iters: dict[str, list[int]] = {m: [] for m in methods}
    secs: dict[str, list[float]] = {m: [] for m in methods}
    disc_err = []
    hist: dict[str, list[float]] = {}
    for n in sizes:
        f, u = manufactured_problem(n)
        for m in methods:
            runs = [solve_poisson(f, m, tol=1e-10) for _ in range(3)]  # best of 3 for timing
            r = min(runs, key=lambda q: q.setup_time + q.solve_time)
            iters[m].append(r.iterations)
            secs[m].append(r.setup_time + r.solve_time)
            if n == sizes[-1] or n == 127:
                hist[(m, n)] = r.residuals  # type: ignore[index]
            if m == "direct":
                disc_err.append(l2_error(r.x, u))
    data["poisson"] = {
        "sizes": sizes,
        "iterations": iters,
        "seconds": secs,
        "discretisation_error": disc_err,
    }

    rows = []
    for m in methods:
        if m == "direct":
            cells = [f"{secs[m][i]:.3f} s (direct)" for i in range(len(sizes))]
        else:
            cells = [f"{iters[m][i]} it / {secs[m][i]:.3f} s" for i in range(len(sizes))]
        rows.append([m, *cells])
    poisson_table = md_table(["Solver", *[f"n = {n}" for n in sizes]], rows)
    hs = np.array([1.0 / (n + 1) for n in sizes])
    d_orders = np.log(np.array(disc_err[:-1]) / np.array(disc_err[1:])) / np.log(hs[:-1] / hs[1:])
    disc_table = md_table(
        ["n", "h", "L2 error of the 5-point solution", "observed order"],
        [
            [str(n), sci(h_, 2), sci(e), "-" if i == 0 else num(float(d_orders[i - 1]))]
            for i, (n, h_, e) in enumerate(zip(sizes, hs, disc_err, strict=True))
        ],
    )

    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    n_hist = 127 if 127 in sizes else sizes[-1]
    for m in ("cg", "cg-ic0", "mg", "mg-cg"):
        res = hist[(m, n_hist)]  # type: ignore[index]
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

    # Heat equation: temporal order, scheme x linear solver
    n = 31
    t_final = 0.1
    ks = [10, 20, 40, 80]
    rows = []
    fig, axes = pl.new_figure(ncols=2, width=5.8, height=4.2)
    heat_data = {}
    for scheme in ("backward-euler", "crank-nicolson"):
        u0, ex = heat_exact_semidiscrete(n, (1, 2), t_final)
        errs = []
        for k in ks:
            hr = heat_solve(u0, t_final, k, scheme, "direct")
            errs.append(float(np.linalg.norm(hr.u - ex) / (n + 1)))
        dts = [t_final / k for k in ks]
        axes[0].loglog(
            dts,
            errs,
            label=scheme,
            **pl.style(scheme if scheme != "crank-nicolson" else "crank-nicolson"),
        )
        o = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(
            np.array(dts[:-1]) / np.array(dts[1:])
        )
        heat_data[scheme] = {"dt": dts, "error": errs, "orders": o.tolist()}
        rows.append([scheme, *[sci(e) for e in errs], num(float(o[-1]))])
    axes[0].set_xlabel("time step dt")
    axes[0].set_ylabel("L2 error at t = 0.1")
    axes[0].set_title("Heat equation: temporal convergence")
    dd = np.array([t_final / 80, t_final / 10])
    pl.loglog_slope_guide(axes[0], dd, 2e-3, 1, "dt", anchor="max")
    pl.loglog_slope_guide(axes[0], dd, 1e-4, 2, "dt^2", anchor="max")
    axes[0].legend()
    heat_table = md_table(
        ["Scheme", *[f"dt = {t_final / k:.4f}" for k in ks], "last-pair order"], rows
    )
    # linear-solver comparison inside the CN time loop
    n2 = 63 if not quick else 31
    u_init = np.zeros((n2, n2))
    u_init[n2 // 3 : n2 // 3 + 6, n2 // 2 : n2 // 2 + 9] = 1.0
    inner_rows = []
    sols = {}
    for ls in ("direct", "cg", "cg-ic0", "mg-cg"):
        hr = min(
            (heat_solve(u_init, 0.02, 20, "crank-nicolson", ls, tol=1e-12) for _ in range(3)),
            key=lambda q: q.wall_time,
        )
        sols[ls] = hr.u
        inner_rows.append(
            [
                ls,
                "1 solve per step" if ls == "direct" else str(hr.total_inner_iterations),
                f"{hr.wall_time:.3f}",
                sci(float(np.max(np.abs(hr.u - sols["direct"])))),
            ]
        )
    heat_solver_table = md_table(
        [
            "Linear solver",
            "inner iterations (20 steps)",
            "wall time (s)",
            "max difference to direct",
        ],
        inner_rows,
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
    data["heat"] = {"orders": heat_data, "solver_comparison": inner_rows}
    return data, {
        "poisson_iterations": poisson_table,
        "poisson_discretisation": disc_table,
        "heat_orders": heat_table,
        "heat_solvers": heat_solver_table,
    }


# --------------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------------


def generate_report(
    out_dir: Path, readme: Path | None = None, quick: bool = False, verbose: bool = True
) -> dict[str, Any]:
    """Run all studies; write ``figures/``, ``results.json``, ``tables.md`` under ``out_dir``."""
    out_dir = Path(out_dir)
    fig_dir = out_dir / "figures"
    res_dir = out_dir / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {"environment": environment(), "quick": quick}
    tables: dict[str, str] = {}
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
            data, tab = fn(fig_dir, quick)
        results[name] = data
        tables.update(tab)
        if verbose:
            print(f"[report] {name}: {time.perf_counter() - t0:.1f} s")
    env = results["environment"]
    tables["environment"] = (
        f"Python {env['python']}, NumPy {env['numpy']}, SciPy {env['scipy']}, "
        f"{env['cpu_count']} CPU cores, {env['platform']}"
    )
    (res_dir / "results.json").write_text(
        json.dumps(results, indent=1, default=float), encoding="utf-8"
    )
    md = "\n\n".join(f"### {k}\n\n{v}" for k, v in tables.items())
    (res_dir / "tables.md").write_text(md + "\n", encoding="utf-8")
    if readme is not None:
        found = update_readme(Path(readme), tables)
        if verbose:
            print(f"[report] README sections updated: {', '.join(found)}")
    return results
