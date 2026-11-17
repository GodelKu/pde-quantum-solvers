"""README tables, docs/results/tables.md and docs/results/results.json must agree.

Every block between ``<!-- BEGIN:name -->`` markers in README.md is rendered from the JSON by
``pdeqs.tables``. Regenerate everything with ``pdeqs report --out docs --readme README.md``, or
refresh the README alone with ``python tools/fill_readme.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from pdeqs.elliptic import heat_exact_semidiscrete, heat_solve, manufactured_problem, solve_poisson
from pdeqs.potentials import DoubleWell, Harmonic, InfiniteWell
from pdeqs.stationary import SOLVERS
from pdeqs.tables import load_results, readme_blocks, render_tables, tables_markdown
from pdeqs.tunneling import transmission_fd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "docs" / "results" / "results.json"


@pytest.fixture(scope="module")
def results() -> dict[str, Any]:
    if not RESULTS.exists():
        pytest.skip("results not generated yet")
    return load_results(RESULTS)


def test_readme_blocks_are_rendered_from_the_results(results):
    blocks = readme_blocks((ROOT / "README.md").read_text(encoding="utf-8"))
    tables = render_tables(results)
    assert set(blocks) == set(tables)
    for name, body in tables.items():
        assert blocks[name] == body, f"README block {name!r} is stale: run tools/fill_readme.py"


def test_tables_md_is_rendered_from_the_results(results):
    text = (ROOT / "docs" / "results" / "tables.md").read_text(encoding="utf-8")
    assert text == tables_markdown(render_tables(results))


def test_results_come_from_a_full_run_and_record_no_host_details(results):
    assert results["quick"] is False
    assert set(results["environment"]) == {"python", "numpy", "scipy", "os"}


# A few deterministic numbers recomputed from the code, so that a change in the numerics
# cannot leave the committed results (and hence the README) silently stale.


def test_stationary_errors_reproduce(results):
    row = results["stationary"]["studies"]["harmonic"]["fd2"]["rows"][0]
    pot = Harmonic()
    e = SOLVERS["fd2"].solve(pot, *pot.box, row["size"], 4).energies
    assert np.max(np.abs(e - pot.exact_energies(4))) == pytest.approx(row["error"], rel=1e-6)

    row = results["stationary"]["hard_walls"]["infinite-well"]["fd4"]["rows"][-1]
    well = InfiniteWell()
    e = SOLVERS["fd4"].solve(well, *well.box, row["size"], 4).energies
    assert np.max(np.abs(e - well.exact_energies(4))) == pytest.approx(row["error"], rel=1e-4)

    dw = results["stationary"]["double_well_splitting"]["chebyshev"]
    e = SOLVERS["chebyshev"].solve(DoubleWell(), *DoubleWell().box, dw["n"], 2).energies
    assert e[1] - e[0] == pytest.approx(dw["splitting"], rel=1e-7)


def test_tunnelling_transmission_reproduces(results):
    tun = results["tunneling"]
    row = next(r for r in tun["fd_convergence"] if r["m"] == 41)
    t = transmission_fd(tun["energy"], 1.0, row["m"], row["h"]).transmission[0]
    assert t == pytest.approx(row["T_fd"], rel=1e-10)


def test_poisson_iteration_counts_reproduce(results):
    po = results["elliptic"]["poisson"]
    n = po["sizes"][0]
    f, _ = manufactured_problem(n)
    for method in ("cg", "cg-ic0", "mg-cg"):
        its = solve_poisson(f, method, tol=po["tolerance"]).iterations
        assert abs(its - po["iterations"][method][0]) <= 1, method  # +-1 for platform rounding


def test_heat_temporal_error_reproduces(results):
    he = results["elliptic"]["heat"]
    n, t_final, k = he["n"], he["t_final"], he["n_steps"][0]
    u0, exact = heat_exact_semidiscrete(n, (1, 2), t_final)
    u = heat_solve(u0, t_final, k, "crank-nicolson", "direct").u
    err = float(np.linalg.norm(u - exact) / (n + 1))
    assert err == pytest.approx(he["orders"]["crank-nicolson"]["error"][0], rel=1e-8)


def test_strang_propagator_splitting_reproduces(results):
    from pdeqs.grid import Grid
    from pdeqs.report import _propagator_splitting

    d = results["timedep"]["double_well_dynamics"]
    row = d["rows"][0]
    dw = DoubleWell()
    grid = Grid.line(d["n"], d["length"])
    (x,) = grid.axes()
    states = SOLVERS["fourier"].solve(dw, -d["length"] / 2, d["length"] / 2, d["n"], 2).states
    split = _propagator_splitting(grid, dw(x), row["dt"], states.astype(complex))
    assert split == pytest.approx(row["splitting_propagator"], rel=1e-9)
