from pathlib import Path

import pytest

from pdeqs.cli import main
from pdeqs.report import md_table, update_readme


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
