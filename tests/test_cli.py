from pathlib import Path

import pytest

from schrodinger_bench.cli import main
from schrodinger_bench.report import md_table, update_readme


def test_convergence_cli_prints_observed_orders(capsys, tmp_path):
    rc = main(["convergence", "crank-nicolson", "--quick", "--out", str(tmp_path)])
    out = capsys.readouterr().out
    assert rc == 0
    assert "observed order" in out and "fitted order" in out
    png = tmp_path / "convergence_crank-nicolson.png"
    assert png.exists() and png.stat().st_size > 5000


@pytest.mark.parametrize("target", ["numerov", "split-step", "poisson", "heat", "tunneling"])
def test_convergence_cli_other_targets(capsys, target):
    assert main(["convergence", target, "--quick", "--no-plot"]) == 0
    assert "error" in capsys.readouterr().out


def test_other_commands(capsys):
    assert main(["list"]) == 0
    assert main(["eigen", "--method", "fourier", "-n", "64", "--states", "3"]) == 0
    assert main(["tunneling", "--energy", "0.72"]) == 0
    assert main(["poisson", "-n", "31", "--solver", "mg-cg"]) == 0
    assert main(["evolve", "--steps", "400", "--t-final", "2", "-n", "128"]) == 0
    out = capsys.readouterr().out
    assert "E_0" in out and "T exact" in out and "mg-cg" in out and "split-step" in out


def test_unknown_method_is_rejected():
    with pytest.raises(SystemExit):
        main(["convergence", "no-such-method"])


def test_update_readme_replaces_only_marked_sections(tmp_path: Path):
    readme = tmp_path / "README.md"
    readme.write_text("a\n<!-- BEGIN:t -->\nold\n<!-- END:t -->\nb\n", encoding="utf-8")
    found = update_readme(readme, {"t": md_table(["x"], [["1"]]), "missing": "zzz"})
    text = readme.read_text(encoding="utf-8")
    assert found == ["t"] and "old" not in text and "| 1 |" in text and "zzz" not in text


@pytest.mark.parametrize(
    "argv",
    [
        ["poisson", "-n", "30"],
        ["poisson", "-n", "31", "--tol", "0"],
        ["eigen", "--method", "fourier", "-n", "4", "--states", "6"],
        ["eigen", "-n", "2", "--states", "1"],
        ["evolve", "--steps", "0"],
        ["tunneling", "--energy", "-1"],
        ["tunneling", "--energy", "5000", "--points", "5"],  # outside the lattice band
    ],
)
def test_invalid_values_give_an_argparse_error_not_a_traceback(argv, capsys):
    with pytest.raises(SystemExit) as exc:
        main(argv)
    assert exc.value.code == 2
    assert "error:" in capsys.readouterr().err


def test_eigen_marks_morse_levels_above_the_bound_states(capsys):
    assert main(["eigen", "--potential", "morse", "--states", "7", "-n", "200"]) == 0
    out = capsys.readouterr().out
    assert "not a bound state" in out and "5 bound states" in out


def test_evolve_reports_an_unstable_rk4_step(capsys):
    assert main(["evolve", "--method", "rk4", "--steps", "50", "-n", "128"]) == 0
    assert "unstable time step" in capsys.readouterr().out
