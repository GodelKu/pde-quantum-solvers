"""Markdown tables for the README, rendered from ``docs/results/results.json`` only.

``sbench report`` runs the studies and writes the JSON; :func:`render_tables` turns that JSON into
the blocks between ``<!-- BEGIN:name -->`` and ``<!-- END:name -->`` in README.md.  Because the
tables depend on nothing but the JSON, ``tests/test_readme_sync.py`` can check that README,
``tables.md`` and the JSON agree without rerunning anything, and ``tools/fill_readme.py``
refreshes the README after a manual edit of the prose.
"""

from __future__ import annotations

import json
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any

# --------------------------------------------------------------------------------------------
# Formatting helpers
# --------------------------------------------------------------------------------------------


def _missing(x: float | None) -> bool:
    return x is None or x != x


def sci(x: float | str | None, digits: int = 2) -> str:
    """Scientific notation; ``None``/NaN print as ``-``, and "inf" (from the JSON) as is."""
    if isinstance(x, str):
        return x
    return "-" if _missing(x) else f"{x:.{digits}e}"


def num(x: float | None, digits: int = 2) -> str:
    return "-" if _missing(x) else f"{x:.{digits}f}"


def timing(stats: dict[str, float], scale: float = 1.0, digits: int = 3) -> str:
    """``median (min-max)`` of repeated wall times, scaled (for example by 1e3 for ms)."""
    lo, mid, hi = (scale * stats[k] for k in ("min", "median", "max"))
    return f"{mid:.{digits}f} ({lo:.{digits}f}-{hi:.{digits}f})"


def md_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def _last_order(study: dict[str, Any]) -> float | None:
    return study["rows"][-1]["order"] if len(study["rows"]) > 1 else None


def _spectral(study: dict[str, Any]) -> bool:
    return study["expected_order"] == "spectral"


def _first_last(study: dict[str, Any]) -> tuple[str, str]:
    first, last = study["rows"][0], study["rows"][-1]
    return (
        f"{sci(first['error'])} (N={first['size']})",
        f"{sci(last['error'])} (N={last['size']})",
    )


# --------------------------------------------------------------------------------------------
# 1. Stationary states
# --------------------------------------------------------------------------------------------


def _stationary(d: dict[str, Any], out: dict[str, str]) -> None:
    for pname, methods in d["studies"].items():
        rows = []
        for method, st in methods.items():
            spectral = _spectral(st)
            first, last = _first_last(st)
            rows.append(
                [
                    method,
                    "spectral" if spectral else st["expected_order"],
                    "-" if spectral else num(st["fitted_order"]),
                    "-" if spectral else num(_last_order(st)),
                    first,
                    last,
                    timing(st["time_largest_n"], 1e3, 1),
                ]
            )
        out[f"stationary_{pname.replace('-', '_')}"] = md_table(
            [
                "Method",
                "Nominal order",
                "Fitted order",
                "Last-pair order",
                "Error, smallest N",
                "Error, largest N",
                "Time at largest N (ms)",
            ],
            rows,
        )

    rows = []
    for pname in ("infinite-well", "linear-ramp"):
        for method, st in d["hard_walls"][pname].items():
            first, last = _first_last(st)
            order = "spectral" if _spectral(st) else num(_last_order(st))
            rows.append([pname, method, first, last, order])
    out["hard_walls"] = md_table(
        ["Problem", "Method", "Error, smallest N", "Error, largest N", "Last-pair order"], rows
    )

    rows = []
    for method, r in d["wavefunctions"].items():
        rows.append([method, sci(r["psi0"]), sci(r["psi3"]), sci(r["energy_6_levels"])])
    out["stationary_wavefunctions"] = md_table(
        ["Method", "max error of psi_0", "max error of psi_3", "max energy error (6 levels)"],
        rows,
    )

    rows = [
        [method, str(r["n"]), f"{r['e0']:.8f}", f"{r['splitting']:.6e}"]
        for method, r in d["double_well_splitting"].items()
    ]
    out["double_well_splitting"] = md_table(
        ["Method", "N", "E_0", "E_1 - E_0 (tunnel splitting)"], rows
    )

    mb = d["morse_box"]
    out["morse_box"] = md_table(
        ["Box", *[f"Chebyshev error, N = {n}" for n in mb["sizes"]]],
        [[f"[{r['left']:g}, {r['right']:g}]", *[sci(e) for e in r["errors"]]] for r in mb["rows"]],
    )

    rows = [
        [
            str(r["n"]),
            timing(r["fd2_eigenpairs"], 1e3, 1),
            timing(r["fd4_eigenvalues"], 1e3, 1),
            timing(r["fd4_eigenpairs"], 1e3, 1),
        ]
        for r in d["fd4_cost"]
    ]
    out["fd4_cost"] = md_table(
        [
            "N",
            "fd2, tridiagonal driver, 4 eigenpairs (ms)",
            "fd4, banded driver, 4 eigenvalues (ms)",
            "fd4, banded driver, 4 eigenpairs (ms)",
        ],
        rows,
    )


# --------------------------------------------------------------------------------------------
# 2. Time-dependent problem
# --------------------------------------------------------------------------------------------


def _timedep(d: dict[str, Any], out: dict[str, str]) -> None:
    rows = []
    for m, st1 in d["temporal_1d"].items():
        st2, sp = d["temporal_2d"][m], d["spatial_1d"][m]
        rows.append(
            [
                m,
                st1["expected_order"],
                num(st1["fitted_order"]),
                num(_last_order(st1)),
                num(st2["fitted_order"]),
                num(_last_order(st2)),
                "spectral" if _spectral(sp) else num(sp["fitted_order"]),
                "-" if _spectral(sp) else num(_last_order(sp)),
            ]
        )
    out["timedep_orders"] = md_table(
        [
            "Method",
            "Nominal dt-order",
            "1D fitted",
            "1D last pair",
            "2D fitted",
            "2D last pair",
            "h-order fitted (1D)",
            "h-order last pair",
        ],
        rows,
    )

    for dim in (1, 2):
        rows = [
            [
                r["method"],
                sci(r["dt"], 1),
                sci(r["max_norm_drift"]),
                sci(r["relative_energy_drift"]),
                sci(r["l2_error_vs_analytic"]),
                timing(r["seconds"]),
            ]
            for r in d[f"conservation_{dim}d"]
        ]
        out[f"conservation_{dim}d"] = md_table(
            [
                "Method",
                "dt",
                "max abs(norm(t) - norm(0))",
                "max rel. energy drift",
                "L2 error vs analytic",
                "Time (s)",
            ],
            rows,
        )

    rk = d["rk4_stability"]
    out["rk4_stability"] = md_table(
        [
            "dt",
            "dt x rho(H)",
            "growth per step of the top mode, abs(R(-i dt rho))",
            f"norm after t = {rk['t_final']:g}",
        ],
        [
            [
                sci(r["dt"], 3),
                num(r["dt_times_rho"], 3),
                num(r["amplification"], 3),
                sci(r["final_norm"], 3),
            ]
            for r in rk["runs"]
        ],
    )

    dw = d["double_well_dynamics"]
    out["double_well_dynamics"] = (
        f"Reference E1 - E0 ({dw['reference']}): {dw['splitting_reference']:.9e}\n\n"
        + md_table(
            [
                "dt",
                "E1 - E0 from the fit of P(x < 0, t)",
                "relative difference",
                "order",
                "E1 - E0 of the Strang step (eigenphases)",
                "relative difference",
                "max abs(norm(t) - norm(0))",
            ],
            [
                [
                    f"{r['dt']:g}",
                    f"{r['splitting_fit']:.9e}",
                    sci(r["relative_difference_fit"]),
                    num(r["order"]),
                    f"{r['splitting_propagator']:.9e}",
                    sci(r["relative_difference_propagator"]),
                    sci(r["max_norm_drift"]),
                ]
                for r in dw["rows"]
            ],
        )
    )

    sc = d["scattering_2d"]
    out["scattering_2d"] = md_table(
        ["Quantity", "Value"],
        [
            ["grid", f"{sc['n']} x {sc['n']}, L = {sc['length']}"],
            ["dt, steps", f"{sc['dt']}, {sc['n_steps']}"],
            ["transmitted probability (split-step)", f"{sc['transmitted_split_step']:.5f}"],
            [
                "transmitted probability (Crank-Nicolson)",
                f"{sc['transmitted_crank_nicolson']:.5f}",
            ],
            [
                "relative L2 difference of the two wave functions",
                sci(sc["relative_l2_difference"]),
            ],
            [
                "same difference with the grid spacing doubled (half as many points)",
                sci(sc["relative_l2_difference_half_grid"]),
            ],
            [
                "max abs(norm(t) - norm(0)): split-step / CN",
                f"{sci(sc['norm_drift_split_step'])} / {sci(sc['norm_drift_crank_nicolson'])}",
            ],
            ["wall time, split-step (s)", timing(sc["seconds_split_step"], digits=2)],
            ["wall time, Crank-Nicolson (s)", timing(sc["seconds_crank_nicolson"], digits=2)],
        ],
    )


# --------------------------------------------------------------------------------------------
# 3. Tunnelling
# --------------------------------------------------------------------------------------------


def _tunneling(d: dict[str, Any], out: dict[str, str]) -> None:
    rows = d["fd_convergence"]
    out["tunneling_fd"] = md_table(
        ["Points in barrier", "h", "T (FD)", "T (exact)", "abs error", "order", "abs(R + T - 1)"],
        [
            [
                str(r["m"]),
                sci(r["h"], 2),
                f"{r['T_fd']:.6f}",
                f"{r['T_exact']:.6f}",
                sci(r["error"]),
                num(r["order"]),
                sci(abs(r["flux_defect"]), 1),
            ]
            for r in rows
        ],
    )
    out["tunneling_packet"] = md_table(
        [
            "k0",
            "mean energy k0^2/2",
            "T(mean energy)",
            "T averaged over packet (exact)",
            "T from split-step",
            "split-step minus averaged exact",
        ],
        [
            [
                num(p["k0"], 2),
                num(0.5 * p["k0"] ** 2, 3),
                f"{p['analytic_mean_energy']:.5f}",
                f"{p['analytic_averaged']:.5f}",
                f"{p['numeric']:.5f}",
                f"{p['numeric'] - p['analytic_averaged']:+.1e}",
            ]
            for p in d["packet"]
        ],
    )
    later = d["packet_later"]
    out["tunneling_packet_later"] = md_table(
        ["k0", "box length", "h", *[f"difference at {m} t1" for m in later[0]["multiples"]]],
        [
            [
                num(r["k0"], 2),
                f"{r['length']:g}",
                f"{r['h']:g}",
                *[f"{x:+.1e}" for x in r["difference"]],
            ]
            for r in later
        ],
    )


# --------------------------------------------------------------------------------------------
# 4. Poisson and heat
# --------------------------------------------------------------------------------------------


def _elliptic(d: dict[str, Any], out: dict[str, str]) -> None:
    po = d["poisson"]
    sizes = po["sizes"]
    heads = [f"n = {n}" for n in sizes]
    rows = [[m, *[str(k) for k in its]] for m, its in po["iterations"].items() if m != "direct"]
    out["poisson_iterations"] = md_table(["Solver", *heads], rows)
    rows = [[m, *[timing(s) for s in secs]] for m, secs in po["seconds"].items()]
    medians = {m: [s["median"] for s in secs] for m, secs in po["seconds"].items()}
    fastest = [min(medians, key=lambda m: medians[m][i]) for i in range(len(sizes))]
    out["poisson_times"] = (
        md_table(["Solver", *[f"{h} (s)" for h in heads]], rows)
        + "\n\nFastest median in this run: "
        + ", ".join(f"{h}: {m}" for h, m in zip(heads, fastest, strict=True))
        + "."
    )
    hs = [1.0 / (n + 1) for n in sizes]
    errs = po["discretisation_error"]
    orders = [None, *po["discretisation_orders"]]
    out["poisson_discretisation"] = md_table(
        ["n", "h", "L2 error of the 5-point solution", "observed order"],
        [
            [str(n), sci(h, 2), sci(e), num(o)]
            for n, h, e, o in zip(sizes, hs, errs, orders, strict=True)
        ],
    )

    he = d["heat"]
    t_final = he["t_final"]
    rows = []
    for scheme, r in he["orders"].items():
        rows.append([scheme, *[sci(e) for e in r["error"]], num(r["orders"][-1])])
    out["heat_orders"] = md_table(
        ["Scheme", *[f"dt = {t_final / k:.4f}" for k in he["n_steps"]], "last-pair order"], rows
    )

    hs_ = he["solvers"]
    steps = hs_["n_steps"]
    rows = []
    for r in hs_["rows"]:
        per_step = "1 LU solve" if r["method"] == "direct" else num(r["iterations"] / steps, 1)
        rows.append(
            [
                r["method"],
                per_step,
                timing(r["setup_seconds"]),
                timing(r["seconds"]),
                sci(r["max_difference_to_direct"]),
            ]
        )
    ranking = sorted(hs_["rows"], key=lambda r: r["seconds"]["median"])
    out["heat_solvers"] = (
        md_table(
            [
                "Linear solver",
                "iterations per step",
                "one-time setup (s)",
                f"total, {steps} steps incl. setup (s)",
                "max difference to direct",
            ],
            rows,
        )
        + "\n\nFastest to slowest median in this run: "
        + ", ".join(r["method"] for r in ranking)
        + "."
    )


# --------------------------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------------------------


def render_tables(results: dict[str, Any]) -> dict[str, str]:
    """Every README block, keyed by marker name, from the results JSON alone."""
    out: dict[str, str] = {}
    _stationary(results["stationary"], out)
    _timedep(results["timedep"], out)
    _tunneling(results["tunneling"], out)
    _elliptic(results["elliptic"], out)
    env = results["environment"]
    out["timing_note"] = (
        f"Wall times are the median of {results['timing']['repeats']} runs after one warm-up "
        "run, with the fastest and slowest of those runs in brackets."
    )
    out["environment"] = (
        f"Python {env['python']}, NumPy {env['numpy']}, SciPy {env['scipy']}, {env['os']}"
    )
    return out


def _pattern(name: str) -> re.Pattern[str]:
    return re.compile(
        rf"(<!-- BEGIN:{re.escape(name)} -->)(.*?)(<!-- END:{re.escape(name)} -->)", re.S
    )


def readme_blocks(text: str) -> dict[str, str]:
    """The current contents of every marked block in ``text``."""
    pattern = re.compile(r"<!-- BEGIN:([\w-]+) -->\n(.*?)\n?<!-- END:\1 -->", re.S)
    return {m.group(1): m.group(2) for m in pattern.finditer(text)}


def fill(text: str, tables: dict[str, str]) -> tuple[str, list[str]]:
    """Replace the marked blocks that exist in ``text``; returns the new text and their names."""
    found = []
    for name, body in tables.items():
        pattern = _pattern(name)
        if pattern.search(text):

            def repl(m: re.Match[str], b: str = body) -> str:
                return f"{m.group(1)}\n{b}\n{m.group(3)}"

            text = pattern.sub(repl, text)
            found.append(name)
    return text, found


def tables_markdown(tables: dict[str, str]) -> str:
    """Content of ``docs/results/tables.md``."""
    return "\n\n".join(f"### {k}\n\n{v}" for k, v in tables.items()) + "\n"


def load_results(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def fill_readme(readme: Path, results: Path) -> bool:
    """Refresh the README blocks from the JSON; returns True if the file changed."""
    text = readme.read_text(encoding="utf-8")
    new, _ = fill(text, render_tables(load_results(results)))
    if new != text:
        readme.write_text(new, encoding="utf-8")
    return new != text
