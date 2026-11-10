import numpy as np
import pytest

from pdeqs.convergence import spatial_study, temporal_study
from pdeqs.grid import Grid
from pdeqs.timedep import (
    RK4,
    CrankNicolson,
    SplitStepFourier,
    coherent_state,
    coherent_state_nd,
    exact_semidiscrete,
    gaussian_packet,
)


def harmonic_problem(n=256, length=20.0, dim=1):
    grid = Grid((n,) * dim, (length,) * dim)
    v = sum(0.5 * m**2 for m in grid.mesh())
    x0 = (1.5, -1.0)[:dim]
    p0 = (0.5, 0.3)[:dim]
    omega = (1.0,) * dim
    return grid, np.asarray(v), (x0, p0, omega)


def test_coherent_state_is_normalised_and_solves_the_equation():
    grid = Grid.line(128, 20.0)
    (x,) = grid.axes()
    psi0 = coherent_state(x, 0.0, 1.0, -0.7)
    assert grid.norm(psi0) == pytest.approx(1.0, abs=1e-12)
    solver = SplitStepFourier(grid, 0.5 * x**2)
    exact = exact_semidiscrete(solver, psi0, 3.0)
    assert grid.norm(exact - coherent_state(x, 3.0, 1.0, -0.7)) < 1e-10


def test_gaussian_packet_is_normalised():
    grid = Grid.square(64, 20.0)
    psi = gaussian_packet(grid, (-3.0, 1.0), 1.2, (1.0, 0.0))
    assert grid.norm(psi) == pytest.approx(1.0, abs=1e-14)


@pytest.mark.parametrize("cls", [SplitStepFourier, CrankNicolson])
def test_unitary_methods_conserve_norm_to_1e_minus_10(cls):
    grid, v, (x0, p0, om) = harmonic_problem()
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    res = cls(grid, v).propagate(psi0, 20.0, 4000)
    assert res.max_norm_drift < 1e-10


def test_split_step_norm_conserved_in_2d():
    grid, v, (x0, p0, om) = harmonic_problem(n=48, length=14.0, dim=2)
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    for cls in (SplitStepFourier, CrankNicolson):
        res = cls(grid, v).propagate(psi0, 4.0, 400)
        assert res.max_norm_drift < 1e-10


def test_crank_nicolson_conserves_energy_to_round_off():
    grid, v, (x0, p0, om) = harmonic_problem()
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    res = CrankNicolson(grid, v).propagate(psi0, 10.0, 500)
    assert res.max_energy_drift / res.energies[0] < 1e-10


def test_split_step_energy_error_scales_as_dt_squared():
    grid, v, (x0, p0, om) = harmonic_problem()
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    drift = [SplitStepFourier(grid, v).propagate(psi0, 4.0, k).max_energy_drift for k in (100, 200)]
    assert 3.5 < drift[0] / drift[1] < 4.5


def test_split_step_matches_analytic_solution():
    grid, v, (x0, p0, om) = harmonic_problem()
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    res = SplitStepFourier(grid, v).propagate(psi0, 2.0, 1000)
    assert grid.norm(res.psi - coherent_state_nd(grid, 2.0, x0, p0, om)) < 1e-5


def test_rk4_loses_norm_exactly_unitary_methods_do_not():
    grid, v, (x0, p0, om) = harmonic_problem(n=128, length=16.0)
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    res = RK4(grid, v).propagate(psi0, 10.0, 1000)  # stable step, coarse on purpose
    assert 1e-12 < res.max_norm_drift < 1e-6


def test_rk4_is_unstable_above_the_stability_limit():
    grid, v, (x0, p0, om) = harmonic_problem(n=256, length=20.0)
    psi0 = coherent_state_nd(grid, 0.0, x0, p0, om)
    solver = RK4(grid, v)
    rho = np.max(np.abs(np.linalg.eigvalsh(solver.h_mat.toarray())))
    steps = int(2.0 * rho / 3.5)  # dt * rho = 3.5 > 2 sqrt(2)
    with np.errstate(all="ignore"):
        res = solver.propagate(psi0, 2.0, steps, n_records=4)
    assert not np.isfinite(res.norms[-1]) or res.norms[-1] > 1e3


@pytest.mark.parametrize(
    ("method", "steps", "low", "high"),
    [
        ("crank-nicolson", [20, 40, 80, 160], 1.9, 2.1),
        ("split-step", [20, 40, 80, 160], 1.95, 2.05),
        ("rk4", [64, 128, 256, 512], 3.9, 4.1),
    ],
)
def test_temporal_order_of_convergence(method, steps, low, high):
    st = temporal_study(method, dim=1, steps=steps)
    assert low < st.fitted_order < high
    assert low < st.pairwise_orders[-1] < high


def test_crank_nicolson_second_order_in_time_2d():
    st = temporal_study("crank-nicolson", dim=2, steps=[20, 40, 80, 160])
    assert 1.9 < st.pairwise_orders[-1] < 2.1


def test_crank_nicolson_second_order_in_space():
    st = spatial_study("crank-nicolson", sizes=[32, 64, 128, 256])
    assert 1.9 < st.fitted_order < 2.1


def test_split_step_spatial_convergence_is_spectral():
    st = spatial_study("split-step", sizes=[16, 24, 32], n_steps=4000)
    assert st.errors[-1] < 1e-6 < st.errors[0] * 1e-1
    # faster than any fixed power: the observed local order keeps growing
    assert st.pairwise_orders[0] > 8
    assert st.pairwise_orders[1] > 1.5 * st.pairwise_orders[0]


def test_propagate_rejects_bad_step_count():
    grid, v, _ = harmonic_problem(n=32)
    with pytest.raises(ValueError):
        CrankNicolson(grid, v).propagate(np.ones(32, dtype=complex), 1.0, 0)
    with pytest.raises(ValueError):
        CrankNicolson(grid, v[:5])
