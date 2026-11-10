import numpy as np
import pytest

from pdeqs.elliptic import (
    Multigrid,
    apply_operator,
    conjugate_gradient,
    heat_exact_semidiscrete,
    heat_solve,
    l2_error,
    manufactured_problem,
    poisson_matrix,
    prolong,
    restrict,
    solve_poisson,
)


def test_matrix_free_operator_matches_sparse_matrix():
    n = 15
    rng = np.random.default_rng(1)
    u = rng.standard_normal((n, n))
    for shift in (0.0, 3.5):
        a = poisson_matrix(n, shift)
        assert np.allclose(a @ u.ravel(), apply_operator(u, shift).ravel())


def test_restriction_is_the_scaled_adjoint_of_prolongation():
    rng = np.random.default_rng(2)
    ec = rng.standard_normal((7, 7))
    r = rng.standard_normal((15, 15))
    assert np.sum(prolong(ec) * r) == pytest.approx(4.0 * np.sum(ec * restrict(r)), rel=1e-12)


def test_prolongation_reproduces_linear_functions():
    nc = 7
    s = np.arange(1, nc + 1) / (nc + 1)
    ec = np.outer(s, np.ones(nc))
    fine = prolong(ec)
    n = 2 * nc + 1
    sf = np.arange(1, n + 1) / (n + 1)
    # exact for linear data away from the far boundary (where the padding is the Dirichlet 0)
    assert np.allclose(fine[:-1, 3], sf[:-1], atol=1e-14)


def test_discretisation_error_is_second_order():
    errs = []
    for n in (15, 31, 63, 127):
        f, u = manufactured_problem(n)
        errs.append(l2_error(solve_poisson(f, "direct").x, u))
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2.0)
    assert np.all(np.abs(orders - 2.0) < 0.05)


@pytest.mark.parametrize("method", ["cg", "cg-jacobi", "cg-ic0", "mg", "mg-cg"])
def test_all_solvers_agree_with_direct_solution(method):
    f, _ = manufactured_problem(31)
    direct = solve_poisson(f, "direct").x
    res = solve_poisson(f, method, tol=1e-11)
    assert res.converged
    assert np.max(np.abs(res.x - direct)) < 1e-9 * np.max(np.abs(direct)) * 1e3


def test_cg_iterations_grow_linearly_with_n_and_jacobi_changes_nothing():
    its, its_j = [], []
    for n in (31, 63, 127):
        f, _ = manufactured_problem(n)
        its.append(solve_poisson(f, "cg").iterations)
        its_j.append(solve_poisson(f, "cg-jacobi").iterations)
    assert its == its_j  # constant diagonal: Jacobi scaling is a no-op for CG
    assert 1.7 < its[1] / its[0] < 2.3 and 1.7 < its[2] / its[1] < 2.3


def test_incomplete_cholesky_beats_plain_cg():
    f, _ = manufactured_problem(63)
    assert solve_poisson(f, "cg-ic0").iterations < 0.5 * solve_poisson(f, "cg").iterations


def test_multigrid_contraction_is_grid_independent():
    rates = []
    for n in (31, 63, 127, 255):
        f, _ = manufactured_problem(n)
        res = solve_poisson(f, "mg", tol=1e-10)
        rates.append(res.mean_reduction)
        assert res.iterations <= 16
    assert max(rates) < 0.25
    assert max(rates) - min(rates) < 0.05


def test_multigrid_preconditioned_cg_iteration_count_is_bounded():
    counts = []
    for n in (31, 63, 127, 255):
        f, _ = manufactured_problem(n)
        counts.append(solve_poisson(f, "mg-cg", tol=1e-10).iterations)
    assert max(counts) <= 10 and max(counts) - min(counts) <= 2


def test_multigrid_requires_power_of_two_grids():
    with pytest.raises(ValueError):
        Multigrid(30)
    with pytest.raises(ValueError):
        solve_poisson(np.zeros((7, 7)), "bogus")


def test_cg_solves_small_spd_system():
    rng = np.random.default_rng(3)
    m = rng.standard_normal((20, 20))
    a = m @ m.T + 20 * np.eye(20)
    b = rng.standard_normal(20)
    res = conjugate_gradient(lambda v: a @ v, b, tol=1e-12)
    assert np.allclose(a @ res.x, b, atol=1e-9)
    assert res.iterations <= 20


def test_heat_equation_temporal_orders():
    n, t_final = 31, 0.1
    u0, exact = heat_exact_semidiscrete(n, (1, 2), t_final)
    for scheme, expected in (("backward-euler", 1.0), ("crank-nicolson", 2.0)):
        errs = []
        for k in (20, 40, 80):
            u = heat_solve(u0, t_final, k, scheme, "mg-cg").u
            errs.append(np.linalg.norm(u - exact) / (n + 1))
        order = np.log(errs[1] / errs[2]) / np.log(2.0)
        assert abs(order - expected) < 0.1, scheme


def test_heat_linear_solvers_agree_and_decay_is_monotone():
    n = 31
    u0 = np.zeros((n, n))
    u0[10:16, 12:21] = 1.0
    sols = {
        ls: heat_solve(u0, 0.01, 10, "crank-nicolson", ls, tol=1e-13).u
        for ls in ("direct", "cg", "mg-cg")
    }
    assert np.max(np.abs(sols["cg"] - sols["direct"])) < 1e-10
    assert np.max(np.abs(sols["mg-cg"] - sols["direct"])) < 1e-10
    assert np.linalg.norm(sols["direct"]) < np.linalg.norm(u0)
