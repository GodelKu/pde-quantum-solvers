import numpy as np
import pytest
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve_triangular

from schrodinger_bench.elliptic import (
    Multigrid,
    PoissonSolver,
    apply_operator,
    conjugate_gradient,
    heat_exact_semidiscrete,
    heat_solve,
    ic0_factor,
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
    assert np.max(np.abs(res.x - direct)) < 1e-6 * np.max(np.abs(direct))


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


def test_heat_solve_sets_up_the_linear_solver_once(monkeypatch):
    """The step matrix is constant: one LU (or one pair of IC(0) triangular factors) per run."""
    import schrodinger_bench.elliptic as el

    calls: list[int] = []
    real_splu = el.splu

    def counting_splu(*args, **kwargs):
        calls.append(1)
        return real_splu(*args, **kwargs)

    monkeypatch.setattr(el, "splu", counting_splu)
    u0, _ = heat_exact_semidiscrete(15, (1, 1), 0.0)
    heat_solve(u0, 0.01, 10, "crank-nicolson", "direct")
    assert len(calls) == 1
    calls.clear()
    heat_solve(u0, 0.01, 10, "crank-nicolson", "cg-ic0")
    assert len(calls) == 2  # L and L^T, factored once each


def test_prepared_solver_matches_one_shot_solve():
    f, _ = manufactured_problem(31)
    for method in ("direct", "cg-ic0", "mg-cg"):
        solver = PoissonSolver(31, method, shift=50.0)
        a = solver.solve(f, tol=1e-12).x
        b = solve_poisson(f, method, shift=50.0, tol=1e-12).x
        assert np.array_equal(a, b), method


@pytest.mark.parametrize("shift", [0.0, 2000.0])
def test_ic0_factor_reproduces_the_matrix_on_its_sparsity_pattern(shift):
    """IC(0): L has the pattern of tril(A), and L L^T = A wherever A is non-zero."""
    n = 8
    a = poisson_matrix(n, shift).toarray()
    low = ic0_factor(n, shift).toarray()
    assert np.all(low[np.tril(a) == 0] == 0.0)
    pattern = a != 0
    assert np.max(np.abs((low @ low.T - a)[pattern])) < 1e-12 * np.max(np.abs(a))
    # the dropped fill-in is what makes it incomplete: L L^T != A off the pattern
    assert np.max(np.abs((low @ low.T)[~pattern])) > 1e-3 * np.max(np.abs(a))


def test_ic0_preconditioner_is_plain_triangular_substitution():
    n, shift = 15, 50.0
    low = sp.csr_matrix(ic0_factor(n, shift))
    r = np.random.default_rng(4).standard_normal(n * n)
    expected = spsolve_triangular(low.T.tocsr(), spsolve_triangular(low, r), lower=False)
    solver = PoissonSolver(n, "cg-ic0", shift)
    assert solver._precond is not None
    assert np.allclose(solver._precond(r), expected, rtol=1e-12, atol=0.0)


def test_iteration_counts_are_pinned():
    """Regression values (README table); +-1 allows for platform rounding at the threshold."""
    expected = {"cg": (106, 215, 434), "cg-ic0": (35, 67, 132), "mg": (13, 13, 14)}
    for i, n in enumerate((31, 63, 127)):
        f, _ = manufactured_problem(n)
        for method, counts in expected.items():
            its = solve_poisson(f, method, tol=1e-10).iterations
            assert abs(its - counts[i]) <= 1, (method, n, its)


@pytest.mark.parametrize("method", ["direct", "cg", "cg-ic0", "mg", "mg-cg"])
def test_zero_right_hand_side_gives_zero_solution(method):
    res = solve_poisson(np.zeros((15, 15)), method)
    assert res.converged and not np.any(res.x)


def test_heat_solve_from_zero_stays_zero():
    for ls in ("direct", "cg", "mg-cg"):
        u = heat_solve(np.zeros((15, 15)), 0.1, 3, "crank-nicolson", ls).u
        assert np.all(u == 0.0), ls


def test_crank_nicolson_undershoots_on_discontinuous_data_backward_euler_does_not():
    """CN is A-stable but not L-stable: the stiff modes of a discontinuous patch flip sign."""
    n = 63
    u0 = np.zeros((n, n))
    u0[n // 3 : n // 3 + 6, n // 2 : n // 2 + 9] = 1.0
    cn = heat_solve(u0, 0.001, 1, "crank-nicolson", "direct").u
    be = heat_solve(u0, 0.001, 1, "backward-euler", "direct").u
    assert cn.min() < -0.05
    assert be.min() > -1e-14  # (I - dt Laplacian_h)^-1 has non-negative entries
