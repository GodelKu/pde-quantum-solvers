"""Stationary 1D Schrodinger equation ``-psi''/2 + V psi = E psi`` (``hbar = m = 1``).

Every solver implements :class:`StationarySolver`, so they can be swapped freely::

    result = SOLVERS["numerov"].solve(potential, x_min, x_max, n_points=400, n_states=6)

``n_points`` is the resolution parameter: the number of interior nodes (finite differences,
Numerov, shooting, Chebyshev) or of Fourier modes (periodic grid).
Dirichlet boundary conditions ``psi(x_min) = psi(x_max) = 0`` are used, except for the
Fourier method, which is periodic.
"""

from __future__ import annotations

import time
import warnings
from abc import ABC, abstractmethod
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, ClassVar, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eig, eig_banded, eigh, eigh_tridiagonal, solve_banded

FloatArray: TypeAlias = NDArray[np.floating[Any]]
Potential = Callable[[FloatArray], FloatArray]


@dataclass
class StationaryResult:
    """Lowest eigenpairs on the solver's own grid ``x`` (states normalised to unit L2 norm)."""

    method: str
    x: FloatArray
    energies: FloatArray
    states: FloatArray  # shape (len(x), n_states)
    n_points: int
    wall_time: float = 0.0
    meta: dict[str, float] = field(default_factory=dict)


class StationarySolver(ABC):
    """Common interface of all stationary-state solvers."""

    name: ClassVar[str]
    order: ClassVar[str]  # nominal convergence order, for documentation

    def solve(
        self,
        potential: Potential,
        x_min: float,
        x_max: float,
        n_points: int,
        n_states: int = 6,
    ) -> StationaryResult:
        if not x_max > x_min:
            raise ValueError("x_max must exceed x_min")
        if n_states < 1 or n_states > n_points:
            raise ValueError("need 1 <= n_states <= n_points")
        t0 = time.perf_counter()
        x, e, psi = self._solve(potential, float(x_min), float(x_max), int(n_points), n_states)
        dt = time.perf_counter() - t0
        return StationaryResult(self.name, x, e, psi, n_points, dt)

    @abstractmethod
    def _solve(
        self, potential: Potential, x_min: float, x_max: float, n: int, k: int
    ) -> tuple[FloatArray, FloatArray, FloatArray]: ...


def _interior_grid(x_min: float, x_max: float, n: int) -> tuple[FloatArray, float]:
    h = (x_max - x_min) / (n + 1)
    return x_min + h * np.arange(1, n + 1), h


def _normalise_columns(psi: FloatArray, weights: FloatArray | float) -> FloatArray:
    norms = np.sqrt(np.sum(psi**2 * np.reshape(weights, (-1, 1)), axis=0))
    return psi / norms


class FiniteDifference2(StationarySolver):
    """Three-point central differences, ``O(h^2)``; tridiagonal symmetric eigenproblem."""

    name = "fd2"
    order = "2"

    def _solve(self, potential, x_min, x_max, n, k):
        x, h = _interior_grid(x_min, x_max, n)
        d = 1.0 / h**2 + potential(x)
        e = np.full(n - 1, -0.5 / h**2)
        w, v = eigh_tridiagonal(d, e, select="i", select_range=(0, k - 1))
        return x, w, _normalise_columns(v, h)


class FiniteDifference4(StationarySolver):
    """Five-point central differences, ``O(h^4)``; symmetric banded eigenproblem."""

    name = "fd4"
    order = "4"

    def _solve(self, potential, x_min, x_max, n, k):
        x, h = _interior_grid(x_min, x_max, n)
        # -1/2 * (-1/12, 4/3, -5/2, 4/3, -1/12)/h^2 : diagonal 5/4, off-diagonals -2/3, +1/24
        ab = np.zeros((3, n))  # lower banded storage
        ab[0] = 1.25 / h**2 + potential(x)
        ab[1, : n - 1] = -(2.0 / 3.0) / h**2
        ab[2, : n - 2] = (1.0 / 24.0) / h**2
        w, v = eig_banded(ab, lower=True, select="i", select_range=(0, k - 1))
        return x, w, _normalise_columns(v, h)


class Numerov(StationarySolver):
    """Matrix Numerov method, ``O(h^4)``.

    Numerov's formula for ``psi'' = 2 (V - E) psi`` reads ``D2 psi = 2 M (V - E) psi`` with
    ``D2`` the three-point second difference and ``M = I + h^2 D2 / 12``.  Since ``M`` and
    ``D2`` commute, ``M^{-1} (-D2/2 + M V) = -M^{-1} D2 / 2 + V`` is a *symmetric* matrix whose
    eigenvalues are exactly the Numerov eigenvalues.
    """

    name = "numerov"
    order = "4"

    def _solve(self, potential, x_min, x_max, n, k):
        x, h = _interior_grid(x_min, x_max, n)
        d2 = (
            np.diag(np.full(n - 1, 1.0), -1) + np.diag(np.full(n - 1, 1.0), 1) - 2.0 * np.eye(n)
        ) / h**2
        m_band = np.zeros((3, n))
        m_band[0, 1:] = 1.0 / 12.0  # M = I + h^2 D2 / 12: off-diagonals 1/12, diagonal 5/6
        m_band[1] = 5.0 / 6.0
        m_band[2, :-1] = 1.0 / 12.0
        minv_d2 = solve_banded((1, 1), m_band, d2)
        kin = -0.5 * 0.5 * (minv_d2 + minv_d2.T)
        ham = kin + np.diag(potential(x))
        w, v = eigh(ham, subset_by_index=(0, k - 1))
        return x, w, _normalise_columns(v, h)


class FourierSpectral(StationarySolver):
    """Fourier pseudospectral collocation on a periodic grid (spectral accuracy).

    The domain is treated as periodic, which is exact for states that vanish at the box edge.
    ``n_points`` should be even.
    """

    name = "fourier"
    order = "spectral"

    def _solve(self, potential, x_min, x_max, n, k):
        length = x_max - x_min
        h = length / n
        x = x_min + h * np.arange(n)
        kk = 2.0 * np.pi * np.fft.fftfreq(n, d=h)
        eye = np.eye(n)
        kin = np.fft.ifft(0.5 * kk[:, None] ** 2 * np.fft.fft(eye, axis=0), axis=0).real
        kin = 0.5 * (kin + kin.T)
        w, v = eigh(kin + np.diag(potential(x)), subset_by_index=(0, k - 1))
        return x, w, _normalise_columns(v, h)


def chebyshev_matrix(n: int) -> tuple[FloatArray, FloatArray]:
    """Chebyshev differentiation matrix on ``n + 1`` Gauss-Lobatto points of ``[-1, 1]``."""
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2.0, np.ones(n - 1), 2.0]) * (-1.0) ** np.arange(n + 1)
    xx = np.tile(x, (n + 1, 1)).T
    dx = xx - xx.T
    d = np.outer(c, 1.0 / c) / (dx + np.eye(n + 1))
    d -= np.diag(d.sum(axis=1))
    return d, x


def clenshaw_curtis_weights(n: int) -> FloatArray:
    """Clenshaw-Curtis quadrature weights for the ``n + 1`` Chebyshev points on ``[-1, 1]``."""
    theta = np.pi * np.arange(n + 1) / n
    w = np.zeros(n + 1)
    ii = np.arange(1, n)
    v = np.ones(n - 1)
    if n % 2 == 0:
        w[0] = w[n] = 1.0 / (n**2 - 1)
        for j in range(1, n // 2):
            v -= 2.0 * np.cos(2 * j * theta[ii]) / (4 * j**2 - 1)
        v -= np.cos(n * theta[ii]) / (n**2 - 1)
    else:
        w[0] = w[n] = 1.0 / n**2
        for j in range(1, (n - 1) // 2 + 1):
            v -= 2.0 * np.cos(2 * j * theta[ii]) / (4 * j**2 - 1)
    w[ii] = 2.0 * v / n
    return w


class ChebyshevSpectral(StationarySolver):
    """Chebyshev collocation with Dirichlet ends (spectral accuracy); dense non-symmetric."""

    name = "chebyshev"
    order = "spectral"

    def _solve(self, potential, x_min, x_max, n, k):
        big_n = n + 1  # n interior nodes -> n + 2 Gauss-Lobatto points
        d, xc = chebyshev_matrix(big_n)
        half = 0.5 * (x_max - x_min)
        mid = 0.5 * (x_max + x_min)
        d2 = (d @ d) / half**2
        y = mid + half * xc
        wq = clenshaw_curtis_weights(big_n) * half
        inner = slice(1, big_n)
        ham = -0.5 * d2[inner, inner] + np.diag(potential(y[inner]))
        w, v = eig(ham)
        order = np.argsort(w.real)[:k]
        e = w.real[order]
        psi = v[:, order].real
        # Chebyshev nodes run from x_max to x_min; report ascending x
        x = y[inner][::-1]
        psi = _normalise_columns(psi, wq[inner])[::-1]
        return x, e, psi


class Shooting(StationarySolver):
    """Shooting with Numerov integration and bisection on the node count (Sturm oscillation).

    For trial energy ``E`` the Numerov recurrence is integrated from the left wall.  The number
    of sign changes of ``psi`` on the grid equals the number of eigenvalues below ``E``; bisecting
    on that count brackets ``E_n`` to machine precision without any eigenvalue solver.  The
    eigenfunction is built by matching a left and a right shot at the point of largest amplitude.
    Pure Python loops: this is the slowest method here by construction.
    """

    name = "shooting"
    order = "4"

    @staticmethod
    def _shoot(g: list[float], h2: float, start: int, stop: int, step: int) -> list[float]:
        """Numerov recurrence for psi'' = g psi between node indices (inclusive), psi=0 at start."""
        out = [0.0] * len(g)
        c = [1.0 - h2 * gi / 12.0 for gi in g]
        phi_prev = 0.0
        out[start] = 0.0
        psi = 1.0
        out[start + step] = psi
        phi = c[start + step] * psi
        i = start + step
        while i != stop:
            nxt = i + step
            phi_next = 2.0 * phi - phi_prev + h2 * g[i] * psi
            psi = phi_next / c[nxt]
            out[nxt] = psi
            phi_prev, phi = phi, phi_next
            if abs(phi) > 1e150:  # rescale to avoid overflow; sign pattern is unaffected
                s = 1e-150
                phi_prev *= s
                phi *= s
                psi *= s
                for j in range(start, nxt + step, step):
                    out[j] *= s
            i = nxt
        return out

    def _count(self, v: list[float], e: float, h2: float) -> int:
        g = [2.0 * (vi - e) for vi in v]
        psi = self._shoot(g, h2, 0, len(g) - 1, 1)
        count = 0
        prev = 0.0
        for val in psi[1:]:
            if val == 0.0:
                continue
            if prev != 0.0 and (val > 0) != (prev > 0):
                count += 1
            prev = val
        return count

    def _energy(self, v: list[float], h2: float, idx: int, lo: float, hi0: float) -> float:
        hi = hi0
        while self._count(v, hi, h2) < idx + 1:
            hi = lo + 2.0 * (hi - lo)
            if hi - lo > 1e8:
                raise RuntimeError("could not bracket eigenvalue")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if self._count(v, mid, h2) >= idx + 1:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-14 * max(1.0, abs(hi)):
                break
        return 0.5 * (lo + hi)

    def _solve(self, potential, x_min, x_max, n, k):
        x, h = _interior_grid(x_min, x_max, n)
        xf = np.concatenate(([x_min], x, [x_max]))
        v = potential(xf).tolist()
        h2 = h * h
        vmin = min(v)
        if h2 * (max(v) - vmin) > 6.0:
            warnings.warn(
                "shooting: h^2 (V_max - V_min) > 6, the Numerov node count is unreliable near "
                "steep walls; refine the grid or shrink the box",
                RuntimeWarning,
                stacklevel=3,
            )
        energies = np.empty(k)
        lo = vmin
        for idx in range(k):
            energies[idx] = self._energy(v, h2, idx, lo, lo + 1.0)
            lo = energies[idx] - 1e-9 * max(1.0, abs(energies[idx]))  # states are ordered
        states = np.empty((n, k))
        for idx, e in enumerate(energies):
            states[:, idx] = self._eigenfunction(v, e, h2)[1:-1]
        return x, energies, _normalise_columns(states, h)

    def _eigenfunction(self, v: list[float], e: float, h2: float) -> FloatArray:
        m = len(v)
        g = [2.0 * (vi - e) for vi in v]
        left = np.array(self._shoot(g, h2, 0, m - 1, 1))
        right = np.array(self._shoot(g, h2, m - 1, 0, -1))
        allowed = np.nonzero(np.asarray(v) < e)[0]
        turn = int(allowed[-1]) if len(allowed) else m // 2
        turn = min(max(turn, 2), m - 3)
        match = 1 + int(np.argmax(np.abs(left[1 : turn + 1])))
        match = min(max(match, 1), m - 2)
        psi = np.empty(m)
        psi[: match + 1] = left[: match + 1]
        psi[match:] = right[match:] * (left[match] / right[match])
        return psi


SOLVERS: dict[str, StationarySolver] = {
    s.name: s
    for s in (
        FiniteDifference2(),
        FiniteDifference4(),
        Numerov(),
        FourierSpectral(),
        ChebyshevSpectral(),
        Shooting(),
    )
}


def align_sign(psi: FloatArray, reference: FloatArray) -> FloatArray:
    """Flip the sign of ``psi`` so that it has non-negative overlap with ``reference``."""
    return psi if float(np.dot(psi, reference)) >= 0.0 else -psi
