import numpy as np
import pytest

from pdeqs.potentials import RectangularBarrier, transmission_rectangular
from pdeqs.tunneling import (
    averaged_transmission,
    packet_transmission,
    snap_width,
    transmission_fd,
)

V0, A = 1.0, 1.55


def test_analytic_limits():
    t = transmission_rectangular(np.array([1e-4, 0.5, 5.0, 50.0]), V0, A)
    assert np.all((t > 0) & (t <= 1.0))
    assert t[0] < 1e-2 and t[-1] > 0.99


def test_analytic_is_continuous_through_the_barrier_top():
    below, at, above = (
        float(transmission_rectangular(e, V0, A)[0]) for e in (1 - 1e-6, 1.0, 1 + 1e-6)
    )
    assert below == pytest.approx(at, rel=1e-4)
    assert above == pytest.approx(at, rel=1e-4)


def test_thick_barrier_asymptotics():
    e, a = 0.4, 8.0
    kappa = np.sqrt(2 * (V0 - e))
    asym = 16 * e * (V0 - e) / V0**2 * np.exp(-2 * kappa * a)
    assert float(transmission_rectangular(e, V0, a)[0]) == pytest.approx(asym, rel=1e-3)


def test_resonances_above_the_barrier_give_unit_transmission():
    # sin(q a) = 0  =>  T = 1
    q = np.pi / A
    e = V0 + 0.5 * q**2
    assert float(transmission_rectangular(e, V0, A)[0]) == pytest.approx(1.0, abs=1e-12)


def test_fd_transmission_conserves_flux_and_matches_exact():
    m = 81
    h = A / m
    energies = np.array([0.2, 0.72, 1.0, 1.5, 2.5])
    res = transmission_fd(energies, V0, m, h)
    exact = RectangularBarrier(V0, m * h).transmission(energies)
    assert np.max(np.abs(res.reflection + res.transmission - 1.0)) < 1e-12
    assert np.max(np.abs(res.transmission - exact)) < 1e-4


def test_fd_transmission_converges_with_second_order():
    errs, hs = [], []
    for m in (21, 41, 81):
        h = A / m
        t = transmission_fd(0.72, V0, m, h).transmission[0]
        errs.append(abs(t - transmission_rectangular(0.72, V0, m * h)[0]))
        hs.append(h)
    order = np.log(errs[1] / errs[2]) / np.log(hs[1] / hs[2])
    assert 1.9 < order < 2.1


def test_snap_width_gives_odd_number_of_points():
    width, m = snap_width(1.5, 0.05)
    assert m % 2 == 1 and width == pytest.approx(m * 0.05)


def test_wave_packet_transmission_matches_momentum_averaged_exact_result():
    p = packet_transmission(1.2, V0, A)
    assert p.norm_drift < 1e-10
    assert abs(p.numeric - p.analytic_averaged) < 1e-3
    # averaging matters: the mean-energy value is visibly different
    assert abs(p.analytic_mean_energy - p.analytic_averaged) > 5e-3
    assert averaged_transmission(1.2, 3.0, V0, p.effective_width) == pytest.approx(
        p.analytic_averaged
    )
