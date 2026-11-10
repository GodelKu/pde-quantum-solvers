"""Uniform Cartesian grids shared by the time-dependent solvers.

Points are ``x_j = -L/2 + j h`` with ``h = L / n`` (``j = 0 .. n-1``).  The same grid is used
by finite-difference solvers (implicit zero Dirichlet values outside the box) and by Fourier
solvers (periodic), so that methods can be compared on identical data.
All quantities are in atomic-like units with ``hbar = m = 1``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.floating[Any]]
ComplexArray: TypeAlias = NDArray[np.complex128]


@dataclass(frozen=True)
class Grid:
    """A uniform grid in one or more dimensions."""

    n: tuple[int, ...]
    length: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.n) != len(self.length):
            raise ValueError("n and length must have the same number of dimensions")
        if any(k < 4 for k in self.n):
            raise ValueError("need at least 4 points per dimension")

    @classmethod
    def line(cls, n: int, length: float) -> Grid:
        return cls((n,), (length,))

    @classmethod
    def square(cls, n: int, length: float) -> Grid:
        return cls((n, n), (length, length))

    @property
    def dim(self) -> int:
        return len(self.n)

    @property
    def shape(self) -> tuple[int, ...]:
        return self.n

    @property
    def h(self) -> tuple[float, ...]:
        return tuple(L / k for L, k in zip(self.length, self.n, strict=True))

    @property
    def cell_volume(self) -> float:
        return float(np.prod(self.h))

    @property
    def size(self) -> int:
        return int(np.prod(self.n))

    def axes(self) -> list[FloatArray]:
        return [
            -L / 2 + h * np.arange(k) for L, h, k in zip(self.length, self.h, self.n, strict=True)
        ]

    def mesh(self) -> list[FloatArray]:
        """Coordinate arrays of shape ``grid.shape`` (``indexing='ij'``)."""
        return list(np.meshgrid(*self.axes(), indexing="ij"))

    def wavenumbers(self) -> list[FloatArray]:
        """Angular wavenumbers of the discrete Fourier transform along each axis."""
        return [2 * np.pi * np.fft.fftfreq(k, d=h) for k, h in zip(self.n, self.h, strict=True)]

    def norm(self, psi: NDArray[np.complexfloating] | FloatArray) -> float:
        """Discrete L2 norm ``sqrt(sum |psi|^2 dV)``."""
        return float(np.sqrt(np.sum(np.abs(psi) ** 2) * self.cell_volume))

    def normalize(self, psi: ComplexArray) -> ComplexArray:
        return psi / self.norm(psi)
