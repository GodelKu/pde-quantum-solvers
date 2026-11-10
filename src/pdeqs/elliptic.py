"""Poisson and heat equations on the unit square with sparse linear algebra and multigrid.

Model problem: ``-Laplacian(u) = f`` on ``(0,1)^2`` with ``u = 0`` on the boundary, discretised by
the five-point stencil on an ``n x n`` interior grid (``n = 2^k - 1``, ``h = 1/(n+1)``).

Solvers (all written here, using SciPy only for sparse matrices and the direct/ILU baselines):

* ``direct``      - sparse LU (SuperLU) baseline,
* ``cg``          - conjugate gradients,
* ``cg-jacobi``   - CG with the diagonal preconditioner (identical to plain CG for this stencil,
                    because the diagonal is constant; kept to show exactly that),
* ``cg-ic0``      - CG with an incomplete-Cholesky (zero fill-in) preconditioner written here,
* ``mg``          - stand-alone geometric multigrid V-cycles,
* ``mg-cg``       - CG preconditioned by one symmetric V-cycle.

The multigrid cycle uses weighted-Jacobi smoothing, full-weighting restriction and bilinear
prolongation.  It also supports a mass shift, ``(sigma I - Laplacian) u = f``, which is what the
implicit time steps of the heat equation need.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, TypeAlias

import numpy as np
import scipy.sparse as sp
from numpy.typing import NDArray
from scipy.sparse.linalg import splu

FloatArray: TypeAlias = NDArray[np.floating[Any]]

# --------------------------------------------------------------------------------------------
# Discrete operators
# --------------------------------------------------------------------------------------------


def poisson_matrix(n: int, shift: float = 0.0) -> sp.csr_matrix:
    """Sparse ``shift * I - Laplacian_h`` on the ``n x n`` interior grid (row-major ordering)."""
    h = 1.0 / (n + 1)
    t = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n))
    eye = sp.identity(n)
    a = (sp.kron(eye, t) + sp.kron(t, eye)) / h**2 + shift * sp.identity(n * n)
    return a.tocsr()


def apply_operator(u: FloatArray, shift: float = 0.0) -> FloatArray:
    """Matrix-free ``(shift I - Laplacian_h) u`` for ``u`` of shape ``(n, n)``."""
    n = u.shape[0]
    h2 = (1.0 / (n + 1)) ** 2
    out = 4.0 * u
    out[1:, :] -= u[:-1, :]
    out[:-1, :] -= u[1:, :]
    out[:, 1:] -= u[:, :-1]
    out[:, :-1] -= u[:, 1:]
    return out / h2 + shift * u


def manufactured_problem(n: int) -> tuple[FloatArray, FloatArray]:
    """Right-hand side and exact solution ``u = X(x) X(y)``, ``X(s) = s (1 - s) exp(s)``.

    ``-Laplacian(u) = f`` with ``f = -(X'' X + X X'')`` and ``X'' = -(3 s + s^2) exp(s)``.  Unlike a
    single sine mode this is *not* a discrete eigenvector, so iterative solvers must work.
    """
    h = 1.0 / (n + 1)
    s = h * np.arange(1, n + 1)
    xs = (s - s**2) * np.exp(s)
    xpp = -(3.0 * s + s**2) * np.exp(s)
    u = np.outer(xs, xs)
    f = -(np.outer(xpp, xs) + np.outer(xs, xpp))
    return f, u


def ic0_factor(n: int, shift: float = 0.0) -> sp.csc_matrix:
    """Incomplete Cholesky factor with zero fill, ``A ~ L L^T``, for the five-point matrix."""
    n2 = n * n
    off = -((n + 1) ** 2)
    diag_a = 4.0 * (n + 1) ** 2 + shift
    d = np.empty(n2)
    for i in range(n2):
        acc = diag_a
        if i % n:
            acc -= (off / d[i - 1]) ** 2
        if i >= n:
            acc -= (off / d[i - n]) ** 2
        d[i] = np.sqrt(acc)
    sub1 = off / d[:-1]
    sub1[(np.arange(1, n2) % n) == 0] = 0.0  # no coupling across grid rows
    subn = off / d[: n2 - n]
    return sp.diags([d, sub1, subn], [0, -1, -n], format="csc")


# --------------------------------------------------------------------------------------------
# Conjugate gradients
# --------------------------------------------------------------------------------------------


@dataclass
class IterationResult:
    x: FloatArray
    iterations: int
    residuals: list[float]  # relative residual norms, entry 0 is the initial residual
    converged: bool
    setup_time: float = 0.0
    solve_time: float = 0.0
    extra: dict[str, float] = field(default_factory=dict)

    @property
    def mean_reduction(self) -> float:
        """Geometric-mean residual reduction per iteration."""
        if self.iterations == 0 or self.residuals[-1] <= 0.0:
            return 0.0
        return float((self.residuals[-1] / self.residuals[0]) ** (1.0 / self.iterations))


def conjugate_gradient(
    matvec: Callable[[FloatArray], FloatArray],
    b: FloatArray,
    precond: Callable[[FloatArray], FloatArray] | None = None,
    tol: float = 1e-10,
    maxiter: int = 5000,
) -> IterationResult:
    """Preconditioned CG; stops when ``||r|| / ||b|| < tol`` (true residual recursion)."""
    t0 = time.perf_counter()
    x = np.zeros_like(b)
    r = b.copy()
    bnorm = float(np.linalg.norm(b))
    z = precond(r) if precond else r
    p = z.copy()
    rz = float(np.vdot(r, z))
    hist = [float(np.linalg.norm(r)) / bnorm]
    k = 0
    while hist[-1] >= tol and k < maxiter:
        ap = matvec(p)
        alpha = rz / float(np.vdot(p, ap))
        x += alpha * p
        r -= alpha * ap
        z = precond(r) if precond else r
        rz_new = float(np.vdot(r, z))
        p = z + (rz_new / rz) * p
        rz = rz_new
        k += 1
        hist.append(float(np.linalg.norm(r)) / bnorm)
    return IterationResult(x, k, hist, hist[-1] < tol, 0.0, time.perf_counter() - t0)


# --------------------------------------------------------------------------------------------
# Geometric multigrid
# --------------------------------------------------------------------------------------------


def restrict(r: FloatArray) -> FloatArray:
    """Full-weighting restriction from ``n = 2 nc + 1`` to ``nc`` interior points."""
    nc = (r.shape[0] - 1) // 2
    rp = np.pad(r, 1)
    c = slice(2, 2 * nc + 1, 2)
    lo = slice(1, 2 * nc, 2)
    hi = slice(3, 2 * nc + 2, 2)
    return (
        4.0 * rp[c, c]
        + 2.0 * (rp[lo, c] + rp[hi, c] + rp[c, lo] + rp[c, hi])
        + rp[lo, lo]
        + rp[lo, hi]
        + rp[hi, lo]
        + rp[hi, hi]
    ) / 16.0


def prolong(ec: FloatArray) -> FloatArray:
    """Bilinear interpolation from ``nc`` to ``n = 2 nc + 1`` interior points."""
    nc = ec.shape[0]
    n = 2 * nc + 1
    ep = np.pad(ec, 1)
    out = np.zeros((n, n))
    out[1::2, 1::2] = ec
    out[0::2, 1::2] = 0.5 * (ep[0 : nc + 1, 1 : nc + 1] + ep[1 : nc + 2, 1 : nc + 1])
    out[1::2, 0::2] = 0.5 * (ep[1 : nc + 1, 0 : nc + 1] + ep[1 : nc + 1, 1 : nc + 2])
    out[0::2, 0::2] = 0.25 * (
        ep[0 : nc + 1, 0 : nc + 1]
        + ep[1 : nc + 2, 0 : nc + 1]
        + ep[0 : nc + 1, 1 : nc + 2]
        + ep[1 : nc + 2, 1 : nc + 2]
    )
    return out


class Multigrid:
    """V-cycle multigrid for ``(shift I - Laplacian_h) u = f`` on an ``n = 2^k - 1`` grid."""

    def __init__(
        self,
        n: int,
        shift: float = 0.0,
        pre: int = 2,
        post: int = 2,
        omega: float = 0.8,
        coarsest: int = 3,
    ) -> None:
        if (n + 1) & n:  # n + 1 must be a power of two
            raise ValueError("multigrid needs n = 2^k - 1")
        self.n = n
        self.shift = shift
        self.pre, self.post, self.omega = pre, post, omega
        self.levels: list[int] = []
        m = n
        while m > coarsest:
            self.levels.append(m)
            m = (m - 1) // 2
        self.levels.append(m)
        self._coarse_lu = splu(poisson_matrix(m, shift).tocsc())

    def _smooth(self, u: FloatArray, f: FloatArray, sweeps: int) -> FloatArray:
        h2 = (1.0 / (u.shape[0] + 1)) ** 2
        inv_diag = 1.0 / (4.0 / h2 + self.shift)
        for _ in range(sweeps):
            u = u + self.omega * inv_diag * (f - apply_operator(u, self.shift))
        return u

    def vcycle(self, u: FloatArray, f: FloatArray, level: int = 0) -> FloatArray:
        if level == len(self.levels) - 1:
            return self._coarse_lu.solve(f.ravel()).reshape(f.shape)
        u = self._smooth(u, f, self.pre)
        res = f - apply_operator(u, self.shift)
        ec = self.vcycle(np.zeros(((u.shape[0] - 1) // 2,) * 2), restrict(res), level + 1)
        u = u + prolong(ec)
        return self._smooth(u, f, self.post)

    def precondition(self, r: FloatArray) -> FloatArray:
        """One V-cycle from a zero initial guess (symmetric positive definite operator)."""
        return self.vcycle(np.zeros_like(r), r)


# --------------------------------------------------------------------------------------------
# Unified driver
# --------------------------------------------------------------------------------------------

POISSON_SOLVERS = ("direct", "cg", "cg-jacobi", "cg-ic0", "mg", "mg-cg")


def solve_poisson(
    f: FloatArray,
    method: str = "mg-cg",
    shift: float = 0.0,
    tol: float = 1e-10,
    maxiter: int = 5000,
    mg: Multigrid | None = None,
) -> IterationResult:
    """Solve ``(shift I - Laplacian_h) u = f`` with the selected method."""
    n = f.shape[0]
    b = f.ravel().astype(float)
    t0 = time.perf_counter()
    if method == "direct":
        a = poisson_matrix(n, shift).tocsc()
        lu = splu(a)
        setup = time.perf_counter() - t0
        t1 = time.perf_counter()
        x = lu.solve(b)
        res = float(np.linalg.norm(b - a @ x) / np.linalg.norm(b))
        return IterationResult(
            x.reshape(n, n), 1, [1.0, res], True, setup, time.perf_counter() - t1
        )

    def mv(v: FloatArray) -> FloatArray:
        return apply_operator(v.reshape(n, n), shift).ravel()

    if method == "cg":
        out = conjugate_gradient(mv, b, None, tol, maxiter)
    elif method == "cg-jacobi":
        inv_d = 1.0 / (4.0 * (n + 1) ** 2 + shift)
        setup = time.perf_counter() - t0
        out = conjugate_gradient(mv, b, lambda r: inv_d * r, tol, maxiter)
        out.setup_time = setup
    elif method == "cg-ic0":
        lower = ic0_factor(n, shift)
        solve_l = splu(lower, permc_spec="NATURAL", diag_pivot_thresh=0.0).solve
        solve_lt = splu(lower.T.tocsc(), permc_spec="NATURAL", diag_pivot_thresh=0.0).solve
        setup = time.perf_counter() - t0
        out = conjugate_gradient(mv, b, lambda r: solve_lt(solve_l(r)), tol, maxiter)
        out.setup_time = setup
    elif method in ("mg", "mg-cg"):
        mg = mg or Multigrid(n, shift)
        setup = time.perf_counter() - t0
        if method == "mg-cg":

            def pc(r: FloatArray) -> FloatArray:
                return mg.precondition(r.reshape(n, n)).ravel()

            out = conjugate_gradient(mv, b, pc, tol, maxiter)
        else:
            out = _multigrid_iteration(mg, f, tol, maxiter)
        out.setup_time = setup
        out.extra["levels"] = float(len(mg.levels))
    else:
        raise ValueError(f"unknown method {method!r}; choose from {POISSON_SOLVERS}")
    if method in ("cg", "cg-jacobi", "cg-ic0", "mg-cg"):
        out.x = out.x.reshape(n, n)
    return out


def _multigrid_iteration(mg: Multigrid, f: FloatArray, tol: float, maxiter: int) -> IterationResult:
    t0 = time.perf_counter()
    u = np.zeros_like(f)
    bnorm = float(np.linalg.norm(f))
    hist = [1.0]
    k = 0
    while hist[-1] >= tol and k < maxiter:
        u = mg.vcycle(u, f)
        hist.append(float(np.linalg.norm(f - apply_operator(u, mg.shift))) / bnorm)
        k += 1
    return IterationResult(u, k, hist, hist[-1] < tol, 0.0, time.perf_counter() - t0)


def l2_error(u: FloatArray, exact: FloatArray) -> float:
    """Grid L2 norm of ``u - exact`` (scaled by ``h``)."""
    h = 1.0 / (u.shape[0] + 1)
    return float(np.linalg.norm(u - exact) * h)


# --------------------------------------------------------------------------------------------
# Heat equation u_t = Laplacian(u)
# --------------------------------------------------------------------------------------------


@dataclass
class HeatResult:
    u: FloatArray
    steps: int
    dt: float
    total_inner_iterations: int
    wall_time: float


def heat_solve(
    u0: FloatArray,
    t_final: float,
    n_steps: int,
    scheme: str = "crank-nicolson",
    linear_solver: str = "mg-cg",
    tol: float = 1e-13,
) -> HeatResult:
    """March ``u_t = Laplacian(u)`` with homogeneous Dirichlet data.

    ``scheme`` is ``"backward-euler"`` (order 1) or ``"crank-nicolson"`` (order 2); each step
    solves ``(I/(theta dt) - Laplacian_h) u_new = rhs`` with the chosen linear solver.
    """
    n = u0.shape[0]
    dt = t_final / n_steps
    theta = {"backward-euler": 1.0, "crank-nicolson": 0.5}[scheme]
    shift = 1.0 / (theta * dt)
    mg = Multigrid(n, shift) if linear_solver in ("mg", "mg-cg") else None
    t0 = time.perf_counter()
    u = u0.copy()
    inner = 0
    for _ in range(n_steps):
        # (I - theta dt L) u_new = (I + (1-theta) dt L) u, with L = Laplacian_h = -A_0
        rhs = u - (1.0 - theta) * dt * apply_operator(u)
        sol = solve_poisson(rhs * shift, linear_solver, shift=shift, tol=tol, mg=mg)
        u = sol.x
        inner += sol.iterations
    return HeatResult(u, n_steps, dt, inner, time.perf_counter() - t0)


def heat_exact_semidiscrete(
    n: int, modes: tuple[int, int], t: float
) -> tuple[FloatArray, FloatArray]:
    """Initial datum ``sin(a pi x) sin(b pi y)`` (a discrete eigenvector) and its exact decay."""
    h = 1.0 / (n + 1)
    a, b = modes
    lam = (4.0 / h**2) * (np.sin(a * np.pi * h / 2) ** 2 + np.sin(b * np.pi * h / 2) ** 2)
    x = h * np.arange(1, n + 1)
    u0 = np.outer(np.sin(a * np.pi * x), np.sin(b * np.pi * x))
    return u0, np.exp(-lam * t) * u0


__all__ = [
    "POISSON_SOLVERS",
    "HeatResult",
    "IterationResult",
    "Multigrid",
    "apply_operator",
    "conjugate_gradient",
    "heat_exact_semidiscrete",
    "heat_solve",
    "ic0_factor",
    "l2_error",
    "manufactured_problem",
    "poisson_matrix",
    "prolong",
    "restrict",
    "solve_poisson",
]
