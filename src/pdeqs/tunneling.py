"""Tunnelling through a rectangular barrier: numerical transmission versus the exact formula.

Two independent numerical routes are provided:

* :func:`transmission_fd` - stationary scattering on a finite-difference lattice with exact
  transparent (lattice plane-wave) boundary conditions, solved as a tridiagonal complex system;
* :func:`packet_transmission` - a Gaussian wave packet propagated with the split-step Fourier
  method; the transmitted probability is compared with the exact ``T(E)`` averaged over the
  packet's momentum distribution.

The barrier width is snapped to an odd multiple of the grid spacing, so that a grid point sits at
the barrier centre and the edges fall *between* grid points.  The analytic reference always uses
this effective width.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import quad
from scipy.linalg import solve_banded

from .grid import Grid
from .potentials import transmission_rectangular
from .timedep import SplitStepFourier, gaussian_packet

FloatArray: TypeAlias = NDArray[np.floating[Any]]


def snap_width(width: float, h: float) -> tuple[float, int]:
    """Nearest odd multiple of ``h`` to ``width``; returns ``(effective_width, n_points)``."""
    m = max(1, round((width / h - 1.0) / 2.0) * 2 + 1)
    return m * h, m


@dataclass
class FDTransmission:
    energies: FloatArray
    transmission: FloatArray
    reflection: FloatArray
    h: float
    effective_width: float


def transmission_fd(
    energies: FloatArray | float,
    height: float,
    n_barrier: int,
    h: float,
    margin: int = 6,
) -> FDTransmission:
    """Finite-difference scattering ``T(E)``, ``R(E)`` for a barrier of ``n_barrier`` grid points.

    Uses the lattice dispersion ``E = (1 - cos(q h)) / h^2`` in the boundary conditions, so that
    the only error is the discretisation of the equation inside the barrier region.
    """
    if n_barrier % 2 == 0:
        raise ValueError("n_barrier must be odd")
    en = np.atleast_1d(np.asarray(energies, dtype=float))
    n = n_barrier + 2 * margin
    j0 = (n + 1) // 2
    x = (np.arange(1, n + 1) - j0) * h
    v = np.where(np.abs(x) < 0.5 * n_barrier * h, height, 0.0)
    t_out = np.empty_like(en)
    r_out = np.empty_like(en)
    off = -0.5 / h**2
    for i, e in enumerate(en):
        c = 1.0 - e * h**2
        if not -1.0 < c < 1.0:
            raise ValueError("energy outside the lattice band; refine the grid")
        q = np.arccos(c) / h
        ph = np.exp(1j * q * h)
        ab = np.zeros((3, n), dtype=complex)
        ab[0, 1:] = off
        ab[2, :-1] = off
        ab[1] = 1.0 / h**2 + v - e
        ab[1, 0] += off * ph
        ab[1, -1] += off * ph
        rhs = np.zeros(n, dtype=complex)
        x0 = x[0] - h
        rhs[0] = -off * np.exp(1j * q * x0) * (1.0 - np.exp(2j * q * h))
        psi = solve_banded((1, 1), ab, rhs)
        t_amp = psi[-1] * np.exp(-1j * q * x[-1])
        r_amp = (psi[0] - np.exp(1j * q * x[0])) * np.exp(1j * q * x[0])
        t_out[i] = abs(t_amp) ** 2
        r_out[i] = abs(r_amp) ** 2
    return FDTransmission(en, t_out, r_out, h, n_barrier * h)


def averaged_transmission(k0: float, sigma: float, height: float, width: float) -> float:
    """``int |phi(k)|^2 T(k^2/2) dk`` for a Gaussian packet (density std ``sigma`` in x)."""
    sk = 1.0 / (2.0 * sigma)

    def integrand(k: float) -> float:
        p = np.exp(-0.5 * ((k - k0) / sk) ** 2) / (sk * np.sqrt(2 * np.pi))
        return float(p * transmission_rectangular(0.5 * k**2, height, width)[0])

    lo, hi = max(1e-6, k0 - 10 * sk), k0 + 10 * sk
    return float(quad(integrand, lo, hi, epsabs=1e-13, epsrel=1e-13, limit=200)[0])


@dataclass
class PacketTransmission:
    k0: float
    sigma: float
    numeric: float
    analytic_averaged: float
    analytic_mean_energy: float
    norm_drift: float
    effective_width: float


def packet_transmission(
    k0: float,
    height: float = 1.0,
    width: float = 1.55,
    sigma: float = 3.0,
    length: float = 200.0,
    n: int = 4000,
    start: float = -30.0,
    dt: float = 0.02,
) -> PacketTransmission:
    """Transmitted probability ``P(x > 0)`` after a head-on collision (split-step Fourier)."""
    grid = Grid.line(n, length)
    (x,) = grid.axes()
    eff_width, _ = snap_width(width, grid.h[0])
    v = np.where(np.abs(x) < 0.5 * eff_width, height, 0.0)
    psi0 = gaussian_packet(grid, (start,), sigma, (k0,))
    # time for the packet centre to reach the same distance beyond the barrier
    t_final = 2.0 * abs(start) / k0
    n_steps = int(np.ceil(t_final / dt))
    solver = SplitStepFourier(grid, v)
    res = solver.propagate(psi0, t_final, n_steps, n_records=2)
    prob = np.abs(res.psi) ** 2
    numeric = float(np.sum(prob[x > 0.5 * eff_width]) * grid.h[0])
    return PacketTransmission(
        k0,
        sigma,
        numeric,
        averaged_transmission(k0, sigma, height, eff_width),
        float(transmission_rectangular(0.5 * k0**2, height, eff_width)[0]),
        res.max_norm_drift,
        eff_width,
    )
