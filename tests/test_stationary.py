import numpy as np
import pytest

from pdeqs.convergence import stationary_study
from pdeqs.potentials import DoubleWell, Harmonic, Morse
from pdeqs.stationary import SOLVERS, align_sign, chebyshev_matrix, clenshaw_curtis_weights

# (n_points, tolerance on the max error of the 4 lowest levels)
HARMONIC_CASES = {
    "fd2": (800, 5e-4),
    "fd4": (400, 1e-5),
    "numerov": (400, 5e-6),
    "shooting": (400, 5e-6),
    "fourier": (80, 1e-11),
    "chebyshev": (80, 1e-10),
}


@pytest.mark.parametrize("method", list(HARMONIC_CASES))
def test_harmonic_oscillator_eigenvalues(method):
    n, tol = HARMONIC_CASES[method]
    pot = Harmonic()
    res = SOLVERS[method].solve(pot, *pot.box, n, 4)
    assert np.max(np.abs(res.energies - pot.exact_energies(4))) < tol


def test_spectral_methods_reach_machine_precision_for_oscillator():
    pot = Harmonic()
    for method in ("fourier", "chebyshev"):
        res = SOLVERS[method].solve(pot, *pot.box, 100, 6)
        assert np.max(np.abs(res.energies - pot.exact_energies(6))) < 1e-11


MORSE_CASES = {
    "fd2": (1600, 2e-3),
    "fd4": (800, 1e-5),
    "numerov": (800, 5e-6),
    "shooting": (800, 5e-6),
    "fourier": (192, 1e-10),
    "chebyshev": (192, 1e-10),
}


@pytest.mark.parametrize("method", list(MORSE_CASES))
def test_morse_levels_match_analytic_formula(method):
    n, tol = MORSE_CASES[method]
    pot = Morse()
    assert pot.n_bound == 5
    res = SOLVERS[method].solve(pot, *pot.box, n, 4)
    assert np.max(np.abs(res.energies - pot.exact_energies(4))) < tol


def test_double_well_methods_agree_with_spectral_reference():
    pot = DoubleWell()
    ref = SOLVERS["chebyshev"].solve(pot, *pot.box, 200, 6).energies
    other = SOLVERS["fourier"].solve(pot, *pot.box, 200, 6).energies
    assert np.max(np.abs(ref - other)) < 1e-10
    for method, tol in (("numerov", 1e-5), ("fd4", 1e-5), ("shooting", 1e-5), ("fd2", 5e-3)):
        e = SOLVERS[method].solve(pot, *pot.box, 800, 6).energies
        assert np.max(np.abs(e - ref)) < tol, method
    # tunnel splitting of the lowest doublet is tiny compared with the level spacing
    assert 0 < ref[1] - ref[0] < 0.02 < ref[2] - ref[1]


def test_shooting_reproduces_numerov_eigenvalues_on_the_same_grid():
    """Node-count bisection solves the same discrete problem as the Numerov matrix."""
    pot = Harmonic()
    a = SOLVERS["shooting"].solve(pot, *pot.box, 300, 5).energies
    b = SOLVERS["numerov"].solve(pot, *pot.box, 300, 5).energies
    assert np.max(np.abs(a - b)) < 1e-9


@pytest.mark.parametrize("method", ["fd2", "fd4", "numerov", "fourier"])
def test_states_are_orthonormal_on_their_grids(method):
    pot = Harmonic()
    res = SOLVERS[method].solve(pot, *pot.box, 160, 4)
    h = res.x[1] - res.x[0]
    gram = res.states.T @ res.states * h
    assert np.max(np.abs(gram - np.eye(4))) < 1e-8


@pytest.mark.parametrize(
    ("method", "tol"),
    [("fourier", 1e-9), ("chebyshev", 1e-9), ("numerov", 1e-5), ("shooting", 1e-5)],
)
def test_oscillator_eigenfunctions(method, tol):
    pot = Harmonic()
    res = SOLVERS[method].solve(pot, *pot.box, 200, 4)
    for k in range(4):
        exact = pot.exact_state(res.x, k)
        psi = align_sign(res.states[:, k], exact)
        assert np.max(np.abs(psi - exact)) < tol


@pytest.mark.parametrize(
    ("method", "low", "high"), [("fd2", 1.9, 2.15), ("fd4", 3.8, 4.2), ("numerov", 3.8, 4.2)]
)
def test_finite_difference_observed_orders(method, low, high):
    st = stationary_study(method, "harmonic", sizes=[40, 80, 160, 320, 640])
    assert low < st.fitted_order < high
    assert low < st.pairwise_orders[-1] < high


def test_shooting_observed_order_is_four():
    st = stationary_study("shooting", "harmonic", sizes=[100, 200, 400, 800])
    assert 3.8 < st.fitted_order < 4.2


def test_shooting_warns_when_grid_too_coarse_for_steep_walls():
    pot = Morse()
    with pytest.warns(RuntimeWarning, match="unreliable"):
        SOLVERS["shooting"].solve(pot, *pot.box, 40, 2)


def test_invalid_arguments():
    with pytest.raises(ValueError):
        SOLVERS["fd2"].solve(Harmonic(), 1.0, -1.0, 100, 3)
    with pytest.raises(ValueError):
        SOLVERS["fd2"].solve(Harmonic(), -5.0, 5.0, 10, 20)
    with pytest.raises(ValueError):
        Morse().exact_energies(9)


def test_clenshaw_curtis_weights_integrate_smooth_functions():
    for n in (16, 17):
        w = clenshaw_curtis_weights(n)
        x = np.cos(np.pi * np.arange(n + 1) / n)
        assert w.sum() == pytest.approx(2.0, abs=1e-13)
        assert np.sum(w * x**4) == pytest.approx(2.0 / 5.0, abs=1e-13)


def test_chebyshev_differentiation_matrix_is_exact_for_polynomials():
    d, x = chebyshev_matrix(12)
    assert np.max(np.abs(d @ x**5 - 5 * x**4)) < 1e-10
