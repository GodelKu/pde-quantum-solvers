"""Model potentials (``hbar = m = 1``), with closed-form results where they exist."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol, TypeAlias

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.floating[Any]]


class Potential1D(Protocol):
    """What the stationary-state studies need from a potential."""

    @property
    def name(self) -> str: ...

    @property
    def box(self) -> tuple[float, float]: ...

    def __call__(self, x: FloatArray) -> FloatArray: ...

    def exact_energies(self, n: int) -> FloatArray | None: ...


@dataclass(frozen=True)
class Harmonic:
    """``V = omega^2 x^2 / 2`` with ``E_n = omega (n + 1/2)``."""

    omega: float = 1.0
    name: str = "harmonic"
    box: tuple[float, float] = (-8.0, 8.0)

    def __call__(self, x: FloatArray) -> FloatArray:
        return 0.5 * self.omega**2 * x**2

    def exact_energies(self, n: int) -> FloatArray:
        return self.omega * (np.arange(n) + 0.5)

    def exact_state(self, x: FloatArray, n: int) -> FloatArray:
        """Normalised eigenfunction ``psi_n`` (physicists' Hermite polynomials)."""
        from numpy.polynomial.hermite import hermval
        from scipy.special import gammaln

        s = np.sqrt(self.omega) * x
        coeffs = np.zeros(n + 1)
        coeffs[n] = 1.0
        log_norm = 0.25 * np.log(self.omega / np.pi) - 0.5 * (n * np.log(2.0) + gammaln(n + 1))
        return np.exp(log_norm - 0.5 * s**2) * hermval(s, coeffs)


@dataclass(frozen=True)
class DoubleWell:
    """``V = lam (x^2 - a^2)^2``; no closed form, a high-resolution reference is used."""

    a: float = 2.0
    lam: float = 0.25
    name: str = "double-well"
    box: tuple[float, float] = (-6.0, 6.0)

    def __call__(self, x: FloatArray) -> FloatArray:
        return self.lam * (x**2 - self.a**2) ** 2

    @property
    def barrier_height(self) -> float:
        return self.lam * self.a**4

    def exact_energies(self, n: int) -> FloatArray | None:
        return None


@dataclass(frozen=True)
class Morse:
    """``V = D (1 - exp(-alpha (x - x_e)))^2``.

    Bound levels: ``E_n = w0 (n+1/2) - (w0 (n+1/2))^2 / (4 D)``, ``w0 = alpha sqrt(2 D)``.
    """

    depth: float = 12.0
    alpha: float = 1.0
    x_eq: float = 0.0
    name: str = "morse"
    box: tuple[float, float] = (-2.0, 18.0)

    def __call__(self, x: FloatArray) -> FloatArray:
        return self.depth * (1.0 - np.exp(-self.alpha * (x - self.x_eq))) ** 2

    @property
    def n_bound(self) -> int:
        lam = np.sqrt(2.0 * self.depth) / self.alpha
        return int(np.floor(lam - 0.5)) + 1

    def exact_energies(self, n: int) -> FloatArray:
        if n > self.n_bound:
            raise ValueError(f"Morse potential has only {self.n_bound} bound states")
        w0 = self.alpha * np.sqrt(2.0 * self.depth)
        e = w0 * (np.arange(n) + 0.5)
        return e - e**2 / (4.0 * self.depth)


@dataclass(frozen=True)
class RectangularBarrier:
    """``V = height`` for ``|x| < width/2`` and ``0`` elsewhere."""

    height: float = 1.0
    width: float = 1.5

    def __call__(self, x: FloatArray) -> FloatArray:
        return np.where(np.abs(x) < 0.5 * self.width, self.height, 0.0)

    def transmission(self, energy: FloatArray | float) -> FloatArray:
        """Exact transmission coefficient for plane waves of energy ``E`` (``m = hbar = 1``)."""
        return transmission_rectangular(energy, self.height, self.width)


def transmission_rectangular(energy: FloatArray | float, height: float, width: float) -> FloatArray:
    """Analytic ``T(E)`` through a rectangular barrier (``E < V0``, ``E = V0`` and ``E > V0``)."""
    e = np.atleast_1d(np.asarray(energy, dtype=float))
    t = np.empty_like(e)
    v0, a = height, width
    below = e < v0 * (1 - 1e-12)
    above = e > v0 * (1 + 1e-12)
    eb = e[below]
    kappa = np.sqrt(2.0 * (v0 - eb))
    t[below] = 1.0 / (1.0 + v0**2 * np.sinh(kappa * a) ** 2 / (4.0 * eb * (v0 - eb)))
    ea = e[above]
    q = np.sqrt(2.0 * (ea - v0))
    t[above] = 1.0 / (1.0 + v0**2 * np.sin(q * a) ** 2 / (4.0 * ea * (ea - v0)))
    eq = ~(below | above)
    t[eq] = 1.0 / (1.0 + 0.5 * v0 * a**2)  # limit E -> V0 (m = hbar = 1)
    return t
