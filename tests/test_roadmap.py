"""Planned extensions from docs/ROADMAP.md. Each test fails until its feature is implemented.

The marker is ``xfail(strict=True)``: once an implementation makes a test pass, the suite turns
red until the marker is removed, so a finished item cannot go unnoticed.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import minimize_scalar

from schrodinger_bench import elliptic, stationary, timedep, tunneling
from schrodinger_bench.grid import Grid
from schrodinger_bench.potentials import Harmonic, Morse, transmission_rectangular


@pytest.mark.xfail(strict=True, reason="roadmap item 1")
def test_richardson_extrapolated_eigenvalues() -> None:
    pot = Harmonic()
    exact = pot.exact_energies(4)
    errs = {}
    for method in ("fd2", "numerov"):
        errs[method] = [
            np.max(np.abs(stationary.richardson_energies(method, pot, *pot.box, n, 4) - exact))
            for n in (79, 159, 319)  # combined with 2n + 1 points: h halved exactly
        ]
    orders = {m: np.log2(np.array(e[:-1]) / np.array(e[1:])) for m, e in errs.items()}
    assert errs["fd2"][1] < 1e-5  # plain fd2 at 319 points: about 2e-3
    assert np.all(np.abs(orders["fd2"] - 4.0) < 0.2)
    assert errs["numerov"][0] < 1e-6
    assert np.all(np.abs(orders["numerov"] - 6.0) < 0.2)


@pytest.mark.xfail(strict=True, reason="roadmap item 2")
def test_imaginary_time_ground_state() -> None:
    grid = Grid.line(256, 20.0)
    (x,) = grid.axes()
    psi0 = grid.normalize(np.exp(-((x - 0.3) ** 2)).astype(complex))
    energy, psi = timedep.imaginary_time_ground_state(grid, 0.5 * x**2, psi0, 0.01, 2000)
    assert abs(energy - 0.5) < 1e-9
    assert grid.norm(psi) == pytest.approx(1.0, abs=1e-12)
    ground = Harmonic().exact_state(x, 0)
    assert abs(np.sum(np.conj(psi) * ground) * grid.h[0]) > 1.0 - 1e-8

    morse = Morse()  # box [-2, 18] as the grid [-10, 10) shifted by 8
    grid = Grid.line(512, 20.0)
    (x,) = grid.axes()
    psi0 = grid.normalize(np.exp(-((x + 8.0) ** 2)).astype(complex))
    energy, _ = timedep.imaginary_time_ground_state(grid, morse(x + 8.0), psi0, 0.01, 2000)
    assert abs(energy - morse.exact_energies(1)[0]) < 1e-6


@pytest.mark.xfail(strict=True, reason="roadmap item 3")
def test_two_dimensional_oscillator_degeneracies() -> None:
    from schrodinger_bench.stationary2d import eigenvalues_2d

    e = eigenvalues_2d(lambda xx, yy: 0.5 * (xx**2 + yy**2), 16.0, 127, 10)
    exact = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4], dtype=float)
    assert np.max(np.abs(e - exact)) < 0.02
    assert abs(e[1] - e[2]) < 1e-9  # (1, 0) and (0, 1): degenerate by the x <-> y symmetry
    assert np.bincount(np.rint(e).astype(int))[1:].tolist() == [1, 2, 3, 4]

    e = eigenvalues_2d(lambda xx, yy: 0.5 * (xx**2 + 4.0 * yy**2), 16.0, 127, 6)
    exact = np.array([1.5, 2.5, 3.5, 3.5, 4.5, 4.5])  # (n_x + 1/2) + 2 (n_y + 1/2)
    assert np.max(np.abs(e - exact)) < 0.02


@pytest.mark.xfail(strict=True, reason="roadmap item 4")
def test_transfer_matrix_single_and_double_barrier() -> None:
    a = 1.55
    energies = np.linspace(0.05, 3.0, 297)
    energies = energies[np.abs(energies - 1.0) > 1e-6]
    t = tunneling.transfer_matrix_transmission(energies, [-a / 2, a / 2], [0.0, 1.0, 0.0])
    assert np.max(np.abs(t - transmission_rectangular(energies, 1.0, a))) < 1e-12

    # two barriers of width 1 and height 1 around a well of width 3
    edges, heights = [-2.5, -1.5, 1.5, 2.5], [0.0, 1.0, 0.0, 1.0, 0.0]

    def transmission(e: float) -> float:
        return float(tunneling.transfer_matrix_transmission(np.array([e]), edges, heights)[0])

    scan = np.linspace(0.005, 0.995, 4000)
    ts = np.array([transmission(e) for e in scan])
    peaks = [j for j in range(1, len(scan) - 1) if ts[j - 1] < ts[j] > ts[j + 1]]
    assert len(peaks) == 2  # two resonances below the barrier top
    found = []
    for j in peaks:
        res = minimize_scalar(
            lambda e: -transmission(e),
            bounds=(scan[j - 1], scan[j + 1]),
            method="bounded",
            options={"xatol": 1e-12},
        )
        assert -res.fun > 1.0 - 1e-9  # resonant tunnelling: complete transmission
        found.append(res.x)
    assert found == pytest.approx([0.234968, 0.954202], abs=1e-4)
    assert transmission(0.1) < 0.01


@pytest.mark.xfail(strict=True, reason="roadmap item 5")
def test_crank_nicolson_with_rannacher_start() -> None:
    n = 63
    u0 = np.zeros((n, n))
    u0[n // 3 : n // 3 + 6, n // 2 : n // 2 + 9] = 1.0
    u1 = elliptic.heat_solve(u0, 0.001, 1, "crank-nicolson", "direct", startup_steps=1).u
    assert u1.min() > -1e-12  # plain CN: about -0.11 here

    m, t_final = 31, 0.1
    smooth, exact = elliptic.heat_exact_semidiscrete(m, (1, 2), t_final)
    errs = [
        np.linalg.norm(
            elliptic.heat_solve(smooth, t_final, k, "crank-nicolson", "direct", startup_steps=1).u
            - exact
        )
        / (m + 1)
        for k in (20, 40, 80)
    ]
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(np.abs(orders - 2.0) < 0.1)  # the start does not cost the second order
