"""Empirical convergence studies: error versus step size, observed order, run time.

The *observed order* between two refinements is ``log(e_i / e_{i+1}) / log(s_i / s_{i+1})`` with
``s`` the step size (``dt`` or ``h``); the *fitted order* is the least-squares slope of
``log e`` against ``log s`` using only points above a round-off floor.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray

from .grid import Grid
from .potentials import DoubleWell, Harmonic, Morse, Potential1D
from .stationary import SOLVERS, StationarySolver
from .timedep import (
    RK4,
    TIME_SOLVERS,
    TimeDependentSolver,
    coherent_state_nd,
    exact_semidiscrete,
)

FloatArray: TypeAlias = NDArray[np.floating[Any]]

ROUNDOFF_FLOOR = 1e-11


@dataclass
class ConvergenceStudy:
    """Errors of one method for a sequence of refinements."""

    name: str
    parameter: str  # "dt" or "h"
    steps: FloatArray  # step sizes, decreasing
    errors: FloatArray
    reference: str
    expected: str
    sizes: list[int] = field(default_factory=list)  # n_steps or n_points per refinement
    wall_times: list[float] = field(default_factory=list)

    @property
    def pairwise_orders(self) -> FloatArray:
        return observed_orders(self.steps, self.errors)

    @property
    def fitted_order(self) -> float:
        if self.is_spectral:
            return float("nan")  # algebraic order is not meaningful for exponential convergence
        return fit_order(self.steps, self.errors)

    @property
    def is_spectral(self) -> bool:
        return self.expected == "spectral"

    def rows(self) -> list[dict[str, float | int | str]]:
        orders = [float("nan"), *self.pairwise_orders.tolist()]
        out = []
        for i, (s, e) in enumerate(zip(self.steps, self.errors, strict=True)):
            row: dict[str, float | int | str] = {
                self.parameter: float(s),
                "error": float(e),
                "order": orders[i],
            }
            if self.sizes:
                row["size"] = self.sizes[i]
            if self.wall_times:
                row["seconds"] = self.wall_times[i]
            out.append(row)
        return out


def observed_orders(steps: Sequence[float] | FloatArray, errors: Sequence[float] | FloatArray):
    s = np.asarray(steps, dtype=float)
    e = np.asarray(errors, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(e[:-1] / e[1:]) / np.log(s[:-1] / s[1:])


def fit_order(
    steps: Sequence[float] | FloatArray,
    errors: Sequence[float] | FloatArray,
    floor: float = ROUNDOFF_FLOOR,
) -> float:
    """Least-squares slope of ``log e`` vs ``log s`` over points with ``e > floor``."""
    s = np.asarray(steps, dtype=float)
    e = np.asarray(errors, dtype=float)
    keep = e > floor
    if keep.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(s[keep]), np.log(e[keep]), 1)[0])


# --------------------------------------------------------------------------------------------
# Stationary eigenvalue convergence
# --------------------------------------------------------------------------------------------

STATIONARY_SIZES: dict[str, list[int]] = {
    "fd2": [40, 80, 160, 320, 640, 1280],
    "fd4": [40, 80, 160, 320, 640, 1280],
    "numerov": [40, 80, 160, 320, 640, 1280],
    "shooting": [200, 400, 800, 1600],
    "fourier": [16, 24, 32, 48, 64, 96],
    "chebyshev": [16, 24, 32, 48, 64, 96],
}

POTENTIALS: dict[str, Potential1D] = {
    "harmonic": Harmonic(),
    "double-well": DoubleWell(),
    "morse": Morse(),
}


def reference_energies(potential: Potential1D, n_states: int) -> tuple[FloatArray, str]:
    """Analytic energies if known, else a converged Chebyshev spectrum (N = 240)."""
    exact = potential.exact_energies(n_states)
    if exact is not None:
        return np.asarray(exact), "analytic"
    ref = SOLVERS["chebyshev"].solve(potential, *potential.box, 240, n_states)
    return ref.energies, "Chebyshev N=240"


def stationary_study(
    method: str,
    potential_name: str = "harmonic",
    sizes: Sequence[int] | None = None,
    n_states: int = 4,
) -> ConvergenceStudy:
    """Max error over the lowest ``n_states`` levels against a trusted reference."""
    solver: StationarySolver = SOLVERS[method]
    potential = POTENTIALS[potential_name]
    ref, ref_name = reference_energies(potential, n_states)
    sizes = list(sizes) if sizes is not None else STATIONARY_SIZES[method]
    x_min, x_max = potential.box
    errors, hs, times = [], [], []
    for n in sizes:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            res = solver.solve(potential, x_min, x_max, n, n_states)
        errors.append(float(np.max(np.abs(res.energies - ref))))
        hs.append((x_max - x_min) / (n + 1))
        times.append(res.wall_time)
    expected = solver.order if solver.order != "spectral" else "spectral"
    return ConvergenceStudy(
        f"{method} / {potential_name}",
        "h",
        np.array(hs),
        np.array(errors),
        ref_name,
        expected,
        list(sizes),
        times,
    )


# --------------------------------------------------------------------------------------------
# Time-dependent convergence
# --------------------------------------------------------------------------------------------

#: Default refinement sequences (number of steps over the whole run).
TIME_STEPS: dict[str, list[int]] = {
    "crank-nicolson": [10, 20, 40, 80, 160, 320],
    "split-step": [10, 20, 40, 80, 160, 320],
    "rk4": [64, 128, 256, 512, 1024],
}


def _coherent_setup(
    dim: int, n: int, length: float
) -> tuple[Grid, FloatArray, tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    grid = Grid((n,) * dim, (length,) * dim)
    x0: tuple[float, ...]
    p0: tuple[float, ...]
    omega: tuple[float, ...]
    if dim == 1:
        x0, p0, omega = (1.5,), (0.5,), (1.0,)
    else:
        x0, p0, omega = (1.5, -1.0), (0.5, 0.3), (1.0, 1.4)
    mesh = grid.mesh()
    v = sum(0.5 * w**2 * m**2 for w, m in zip(omega, mesh, strict=True))
    return grid, np.asarray(v, dtype=float), x0, p0, omega


def temporal_study(
    method: str,
    dim: int = 1,
    n_points: int | None = None,
    t_final: float = 2.0,
    steps: Sequence[int] | None = None,
) -> ConvergenceStudy:
    """Time-step convergence against the exact propagator of the method's own spatial operator.

    A coarse grid is used on purpose (dense eigendecomposition of ``H``); because the reference
    shares the spatial discretisation, the measured error is purely temporal.
    """
    n = n_points if n_points is not None else (64 if dim == 1 else 20)
    length = 16.0 if dim == 1 else 12.0
    grid, v, x0, p0, omega = _coherent_setup(dim, n, length)
    solver = TIME_SOLVERS[method](grid, v)
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, omega)
    ref = exact_semidiscrete(solver, psi0, t_final)
    steps = list(steps) if steps is not None else TIME_STEPS[method]
    errors, dts, wall = [], [], []
    for k in steps:
        res = solver.propagate(psi0, t_final, k, n_records=2)
        errors.append(grid.norm(res.psi - ref))
        dts.append(t_final / k)
        wall.append(res.wall_time)
    return ConvergenceStudy(
        f"{method} / time, {dim}D",
        "dt",
        np.array(dts),
        np.array(errors),
        "exact exp(-iHt) of the same spatial operator",
        str(TIME_SOLVERS[method].order),
        list(steps),
        wall,
    )


SPACE_SIZES: dict[str, list[int]] = {
    "crank-nicolson": [32, 64, 128, 256, 512],
    "rk4": [32, 64, 128, 256, 512],
    "split-step": [12, 16, 20, 24, 28, 32, 40],
}


#: Time steps used in the spatial studies (small enough not to pollute the spatial error).
SPACE_STEPS: dict[str, int] = {"crank-nicolson": 2000, "rk4": 2000, "split-step": 16000}


def spatial_study(
    method: str,
    sizes: Sequence[int] | None = None,
    t_final: float = 1.0,
    n_steps: int | None = None,
    length: float = 16.0,
) -> ConvergenceStudy:
    """Grid convergence (1D harmonic oscillator) against the analytic coherent state.

    ``dt = t_final / n_steps`` is small enough that the temporal error is far below the spatial
    one over the plotted range.
    """
    sizes = list(sizes) if sizes is not None else SPACE_SIZES[method]
    n_steps = n_steps if n_steps is not None else SPACE_STEPS[method]
    errors, hs, wall = [], [], []
    for n in sizes:
        grid, v, x0, p0, omega = _coherent_setup(1, n, length)
        solver = TIME_SOLVERS[method](grid, v)
        psi0 = coherent_state_nd(grid, 0.0, x0, p0, omega)
        res = solver.propagate(psi0, t_final, n_steps, n_records=2)
        exact = coherent_state_nd(grid, t_final, x0, p0, omega)
        errors.append(grid.norm(res.psi - exact))
        hs.append(grid.h[0])
        wall.append(res.wall_time)
    expected = "spectral" if method == "split-step" else "2"
    return ConvergenceStudy(
        f"{method} / space, 1D",
        "h",
        np.array(hs),
        np.array(errors),
        "analytic coherent state",
        expected,
        list(sizes),
        wall,
    )


# --------------------------------------------------------------------------------------------
# Conservation / cost comparison
# --------------------------------------------------------------------------------------------


def conservation_comparison(
    dim: int = 1, t_final: float = 10.0, n_steps: int = 5000, n_points: int | None = None
) -> list[dict[str, Any]]:
    """Run all propagators with the same ``dt`` on the harmonic coherent-state problem."""
    if dim == 1:
        n, length = (n_points or 512), 30.0
    else:
        n, length = (n_points or 96), 16.0
    grid, v, x0, p0, omega = _coherent_setup(dim, n, length)
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, omega)
    exact = coherent_state_nd(grid, t_final, x0, p0, omega)
    rows: list[dict[str, Any]] = []
    for name, cls in TIME_SOLVERS.items():
        solver: TimeDependentSolver = cls(grid, v)
        res = solver.propagate(psi0, t_final, n_steps)
        rows.append(
            {
                "method": name,
                "dt": res.dt,
                "max_norm_drift": res.max_norm_drift,
                "max_energy_drift": res.max_energy_drift,
                "relative_energy_drift": res.max_energy_drift / abs(res.energies[0]),
                "l2_error_vs_analytic": grid.norm(res.psi - exact),
                "seconds": res.wall_time,
                "setup_seconds": res.setup_time,
            }
        )
    return rows


def rk4_stability_demo(n_points: int = 512, length: float = 30.0, t_final: float = 2.0):
    """Final norm of RK4 for a stable and an unstable time step (limit ``dt < 2.83 / rho(H)``)."""
    grid, v, x0, p0, omega = _coherent_setup(1, n_points, length)
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, omega)
    solver = RK4(grid, v)
    rho = float(np.max(np.abs(np.linalg.eigvalsh(solver.h_mat.toarray()))))
    rows = []
    for k in (int(np.ceil(t_final * rho / 2.0)), int(np.floor(t_final * rho / 3.5))):
        with np.errstate(over="ignore", invalid="ignore"):
            res = solver.propagate(psi0, t_final, k, n_records=4)
        rows.append(
            {
                "dt": t_final / k,
                "dt_times_rho": t_final / k * rho,
                "final_norm": float(res.norms[-1]),
            }
        )
    return rho, rows


def cn_study_quick(dim: int = 1) -> ConvergenceStudy:
    """Small, fast CN temporal study used by the CLI smoke test."""
    return temporal_study("crank-nicolson", dim=dim, steps=[10, 20, 40, 80])
