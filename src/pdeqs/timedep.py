"""Time-dependent Schrodinger equation ``i dpsi/dt = H psi``, ``H = -Laplacian/2 + V`` (1D/2D).

Three propagators share one interface (:class:`TimeDependentSolver`):

* :class:`CrankNicolson` - implicit midpoint rule + 3/5-point finite differences
  (second order in time and space, exactly unitary),
* :class:`SplitStepFourier` - Strang splitting ``e^{-iVdt/2} e^{-iTdt} e^{-iVdt/2}`` with FFTs
  (second order in time, spectral in space, exactly unitary),
* :class:`RK4` - classical Runge-Kutta on the finite-difference Hamiltonian
  (fourth order in time, *not* unitary, conditionally stable).
"""

from __future__ import annotations

import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import ClassVar

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray
from scipy.sparse.linalg import splu

from .grid import ComplexArray, FloatArray, Grid

# --------------------------------------------------------------------------------------------
# Initial states with known time evolution
# --------------------------------------------------------------------------------------------


def coherent_state(
    x: FloatArray, t: float, x0: float = 0.0, p0: float = 0.0, omega: float = 1.0
) -> ComplexArray:
    """Exact solution for a Gaussian coherent state in ``V = omega^2 x^2 / 2``.

    ``psi(x, 0)`` is a unit-width ground-state Gaussian displaced to ``(x0, p0)``; it then follows
    the classical orbit, including the global phase ``exp(-i omega t / 2)``.
    """
    s = np.sqrt(omega) * x
    tau = omega * t
    alpha0 = (np.sqrt(omega) * x0 + 1j * p0 / np.sqrt(omega)) / np.sqrt(2.0)
    alpha = alpha0 * np.exp(-1j * tau)
    amp = np.exp(-0.5 * s**2 + np.sqrt(2.0) * alpha * s - 0.5 * alpha**2 - 0.5 * abs(alpha) ** 2)
    return (omega**0.25 * np.pi**-0.25) * amp * np.exp(-0.5j * tau)


def coherent_state_nd(
    grid: Grid, t: float, x0: tuple[float, ...], p0: tuple[float, ...], omega: tuple[float, ...]
) -> ComplexArray:
    """Product of 1D coherent states (separable harmonic oscillator, any dimension)."""
    psi = np.ones(grid.shape, dtype=complex)
    for ax, x in enumerate(grid.axes()):
        shape = [1] * grid.dim
        shape[ax] = -1
        psi = psi * coherent_state(x, t, x0[ax], p0[ax], omega[ax]).reshape(shape)
    return psi


def gaussian_packet(grid: Grid, center: tuple[float, ...], sigma: float, k0: tuple[float, ...]):
    """Normalised Gaussian packet ``exp(-(x-x0)^2/(4 sigma^2) + i k0 x)`` (density std sigma)."""
    psi = np.ones(grid.shape, dtype=complex)
    for ax, x in enumerate(grid.axes()):
        shape = [1] * grid.dim
        shape[ax] = -1
        f = (2 * np.pi * sigma**2) ** -0.25 * np.exp(
            -((x - center[ax]) ** 2) / (4 * sigma**2) + 1j * k0[ax] * x
        )
        psi = psi * f.reshape(shape)
    return grid.normalize(psi)


# --------------------------------------------------------------------------------------------
# Solvers
# --------------------------------------------------------------------------------------------


@dataclass
class TimeResult:
    """Trajectory summary recorded during :meth:`TimeDependentSolver.propagate`."""

    method: str
    times: FloatArray
    norms: FloatArray
    energies: FloatArray
    psi: ComplexArray
    n_steps: int
    dt: float
    wall_time: float  # seconds spent inside ``step`` (setup/factorisation counted separately)
    setup_time: float

    @property
    def max_norm_drift(self) -> float:
        return float(np.max(np.abs(self.norms - self.norms[0])))

    @property
    def max_energy_drift(self) -> float:
        return float(np.max(np.abs(self.energies - self.energies[0])))


class TimeDependentSolver(ABC):
    """Fixed-step propagator for a time-independent potential given on a :class:`Grid`."""

    name: ClassVar[str]
    order: ClassVar[int]
    unitary: ClassVar[bool]

    def __init__(self, grid: Grid, potential: FloatArray) -> None:
        v = np.asarray(potential, dtype=float)
        if v.shape != grid.shape:
            raise ValueError(f"potential has shape {v.shape}, grid expects {grid.shape}")
        self.grid = grid
        self.v = v

    @abstractmethod
    def hamiltonian(self, psi: ComplexArray) -> ComplexArray:
        """Apply the spatially discretised Hamiltonian used by this method."""

    @abstractmethod
    def prepare(self, dt: float) -> None:
        """Pre-compute everything that depends on the time step (factorisations, phases)."""

    @abstractmethod
    def step(self, psi: ComplexArray) -> ComplexArray:
        """Advance ``psi`` by the prepared ``dt``."""

    def energy(self, psi: ComplexArray) -> float:
        """Expectation value ``<psi|H|psi>`` with the method's own spatial operator."""
        return float(np.real(np.vdot(psi, self.hamiltonian(psi))) * self.grid.cell_volume)

    def dense_hamiltonian(self) -> ComplexArray:
        """Dense matrix of :meth:`hamiltonian` (small grids only; used for exact references)."""
        n = self.grid.size
        mat = np.zeros((n, n), dtype=complex)
        for j in range(n):
            e = np.zeros(n, dtype=complex)
            e[j] = 1.0
            mat[:, j] = self.hamiltonian(e.reshape(self.grid.shape)).ravel()
        return mat

    def propagate(
        self, psi0: ComplexArray, t_final: float, n_steps: int, n_records: int = 100
    ) -> TimeResult:
        """Integrate to ``t_final`` in ``n_steps`` equal steps, recording norm and energy."""
        if n_steps < 1:
            raise ValueError("n_steps must be >= 1")
        dt = t_final / n_steps
        t0 = time.perf_counter()
        self.prepare(dt)
        setup = time.perf_counter() - t0

        every = max(1, n_steps // n_records)
        psi = np.array(psi0, dtype=complex)
        times = [0.0]
        norms = [self.grid.norm(psi)]
        energies = [self.energy(psi)]
        spent = 0.0
        done = 0
        while done < n_steps:
            block = min(every, n_steps - done)
            t1 = time.perf_counter()
            for _ in range(block):
                psi = self.step(psi)
            spent += time.perf_counter() - t1
            done += block
            times.append(done * dt)
            norms.append(self.grid.norm(psi))
            energies.append(self.energy(psi))
        return TimeResult(
            self.name,
            np.array(times),
            np.array(norms),
            np.array(energies),
            psi,
            n_steps,
            dt,
            spent,
            setup,
        )


def fd_laplacian(grid: Grid) -> sp.csr_matrix:
    """Sparse second-difference Laplacian with zero Dirichlet values just outside the grid."""
    total = None
    for ax, (n, h) in enumerate(zip(grid.n, grid.h, strict=True)):
        d = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n)) / h**2
        left = sp.identity(int(np.prod(grid.n[:ax])), format="csr")
        right = sp.identity(int(np.prod(grid.n[ax + 1 :])), format="csr")
        term = sp.kron(sp.kron(left, d), right)
        total = term if total is None else total + term
    assert total is not None
    return total.tocsr()


class _FiniteDifferenceSolver(TimeDependentSolver):
    def __init__(self, grid: Grid, potential: FloatArray) -> None:
        super().__init__(grid, potential)
        self.h_mat: sp.csr_matrix = (
            -0.5 * fd_laplacian(grid) + sp.diags(self.v.ravel(), format="csr")
        ).tocsr()

    def hamiltonian(self, psi: ComplexArray) -> ComplexArray:
        return (self.h_mat @ psi.ravel()).reshape(self.grid.shape)


class CrankNicolson(_FiniteDifferenceSolver):
    """``(1 + i dt H/2) psi_{n+1} = (1 - i dt H/2) psi_n`` with a sparse LU factorisation."""

    name = "crank-nicolson"
    order = 2
    unitary = True

    def prepare(self, dt: float) -> None:
        ident = sp.identity(self.grid.size, format="csc")
        self._rhs = (ident - 0.5j * dt * self.h_mat).tocsr()
        self._lu = splu((ident + 0.5j * dt * self.h_mat).tocsc())

    def step(self, psi: ComplexArray) -> ComplexArray:
        return self._lu.solve(self._rhs @ psi.ravel()).reshape(self.grid.shape)


class RK4(_FiniteDifferenceSolver):
    """Classical fourth-order Runge-Kutta for ``dpsi/dt = -i H psi`` (finite differences).

    Stable only for ``dt * rho(H) < 2 sqrt(2)``; the norm is not conserved exactly.
    """

    name = "rk4"
    order = 4
    unitary = False

    def prepare(self, dt: float) -> None:
        self._dt = dt

    def step(self, psi: ComplexArray) -> ComplexArray:
        y = psi.ravel()
        a = self.h_mat
        dt = self._dt
        k1 = -1j * (a @ y)
        k2 = -1j * (a @ (y + 0.5 * dt * k1))
        k3 = -1j * (a @ (y + 0.5 * dt * k2))
        k4 = -1j * (a @ (y + dt * k3))
        return (y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)).reshape(self.grid.shape)


class SplitStepFourier(TimeDependentSolver):
    """Strang splitting with FFT kinetic propagator (periodic boundary conditions)."""

    name = "split-step"
    order = 2
    unitary = True

    def __init__(self, grid: Grid, potential: FloatArray) -> None:
        super().__init__(grid, potential)
        ks = grid.wavenumbers()
        k2 = np.zeros(grid.shape)
        for ax, k in enumerate(ks):
            shape = [1] * grid.dim
            shape[ax] = -1
            k2 = k2 + (k**2).reshape(shape)
        self._k2 = k2

    def hamiltonian(self, psi: ComplexArray) -> ComplexArray:
        kin = np.fft.ifftn(0.5 * self._k2 * np.fft.fftn(psi))
        return kin + self.v * psi

    def prepare(self, dt: float) -> None:
        self._half_v = np.exp(-0.5j * dt * self.v)
        self._kin = np.exp(-0.5j * dt * self._k2)

    def step(self, psi: ComplexArray) -> ComplexArray:
        psi = self._half_v * psi
        psi = np.fft.ifftn(self._kin * np.fft.fftn(psi))
        return self._half_v * psi


TIME_SOLVERS: dict[str, type[TimeDependentSolver]] = {
    cls.name: cls for cls in (CrankNicolson, SplitStepFourier, RK4)
}


def exact_semidiscrete(
    solver: TimeDependentSolver, psi0: ComplexArray, t: float
) -> NDArray[np.complex128]:
    """Exact ``exp(-i H t) psi0`` for the solver's *own* spatial operator (dense ``eigh``).

    Comparing against this isolates the time-integration error from the spatial error.
    """
    h = solver.dense_hamiltonian()
    h = 0.5 * (h + h.conj().T)
    w, u = np.linalg.eigh(h)
    coeff = u.conj().T @ psi0.ravel()
    return (u @ (np.exp(-1j * w * t) * coeff)).reshape(solver.grid.shape)
