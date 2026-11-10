"""Command-line interface: ``pdeqs <command> ...`` (also ``python -m pdeqs``)."""

from __future__ import annotations

import argparse
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from . import __version__
from .convergence import (
    POTENTIALS,
    STATIONARY_SIZES,
    TIME_STEPS,
    ConvergenceStudy,
    conservation_comparison,
    spatial_study,
    stationary_study,
    temporal_study,
)
from .elliptic import (
    POISSON_SOLVERS,
    heat_exact_semidiscrete,
    heat_solve,
    l2_error,
    manufactured_problem,
    solve_poisson,
)
from .potentials import RectangularBarrier
from .report import format_study, generate_report, num, sci
from .stationary import SOLVERS
from .tunneling import transmission_fd

TIME_METHODS = tuple(TIME_STEPS)
CONVERGENCE_TARGETS = (*STATIONARY_SIZES, *TIME_METHODS, "poisson", "heat", "tunneling")


def _plot_studies(studies: Sequence[ConvergenceStudy], path: Path, title: str) -> Path:
    from . import plotting as pl

    fig, ax = pl.new_figure(width=6.0, height=4.4)
    for k, st in enumerate(studies):
        color = ["#2a78d6", "#eb6834", "#1baf7a", "#6250d6"][k % 4]
        ax.loglog(
            st.steps,
            np.maximum(st.errors, 1e-16),
            marker="os^v"[k % 4],
            color=color,
            lw=1.6,
            label=f"{st.name} ({st.expected})",
        )
    ax.set_xlabel(studies[0].parameter)
    ax.set_ylabel("error")
    ax.set_title(title)
    ax.legend()
    return pl.save(fig, path)


def cmd_convergence(args: argparse.Namespace) -> int:
    studies: list[ConvergenceStudy] = []
    m = args.method
    quick = args.quick
    if m in STATIONARY_SIZES:
        sizes = STATIONARY_SIZES[m][:4] if quick else None
        names = [args.potential] if args.potential != "all" else list(POTENTIALS)
        studies = [stationary_study(m, p, sizes=sizes) for p in names]
    elif m in TIME_METHODS:
        steps = TIME_STEPS[m][:4] if quick else None
        studies.append(temporal_study(m, dim=args.dim, steps=steps))
        if args.dim == 1:
            sizes = [32, 64, 128] if (quick and m != "split-step") else None
            studies.append(
                spatial_study(
                    m, sizes=[12, 16, 20, 24, 28, 32] if quick and m == "split-step" else sizes
                )
            )
    elif m == "poisson":
        sizes = [15, 31, 63] if quick else [15, 31, 63, 127, 255]
        errs = []
        for n in sizes:
            f, u = manufactured_problem(n)
            errs.append(l2_error(solve_poisson(f, "direct").x, u))
        h_arr = np.array([1.0 / (n + 1) for n in sizes])
        studies.append(
            ConvergenceStudy(
                "poisson 5-point / space", "h", h_arr, np.array(errs), "exact solution", "2", sizes
            )
        )
    elif m == "heat":
        n = 31
        ks = [10, 20, 40] if quick else [10, 20, 40, 80, 160]
        for scheme in ("backward-euler", "crank-nicolson"):
            u0, ex = heat_exact_semidiscrete(n, (1, 2), 0.1)
            errs = [
                float(np.linalg.norm(heat_solve(u0, 0.1, k, scheme, "mg-cg").u - ex) / (n + 1))
                for k in ks
            ]
            studies.append(
                ConvergenceStudy(
                    f"heat {scheme} / time",
                    "dt",
                    np.array([0.1 / k for k in ks]),
                    np.array(errs),
                    "exact semi-discrete decay",
                    "1" if scheme == "backward-euler" else "2",
                    ks,
                )
            )
    elif m == "tunneling":
        barrier = RectangularBarrier()
        ms = [5, 11, 21, 41] if quick else [5, 11, 21, 41, 81, 161]
        h_tun = np.array([barrier.width / k for k in ms])
        t_exact = float(barrier.transmission(0.72)[0])
        errs = [
            abs(
                float(transmission_fd(0.72, barrier.height, k, barrier.width / k).transmission[0])
                - t_exact
            )
            for k in ms
        ]
        studies.append(
            ConvergenceStudy(
                "tunneling FD scattering / space",
                "h",
                h_tun,
                np.array(errs),
                "exact T(E=0.72)",
                "2",
                ms,
            )
        )
    for st in studies:
        print(format_study(st))
        print()
    if not args.no_plot and studies:
        out = Path(args.out) / f"convergence_{m}.png"
        print(f"figure written to {_plot_studies(studies, out, f'Convergence: {m}')}")
    return 0


def cmd_eigen(args: argparse.Namespace) -> int:
    pot = POTENTIALS[args.potential]
    x_min, x_max = pot.box
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        res = SOLVERS[args.method].solve(pot, x_min, x_max, args.n_points, args.states)
    exact = (
        pot.exact_energies(args.states) if args.states <= getattr(pot, "n_bound", 10**9) else None
    )
    print(
        f"{args.method} on {args.potential}, N={args.n_points}, box=[{x_min}, {x_max}], "
        f"{res.wall_time * 1e3:.1f} ms"
    )
    for k, e in enumerate(res.energies):
        line = f"  E_{k} = {e:.10f}"
        if exact is not None:
            line += f"   exact {exact[k]:.10f}   error {abs(e - exact[k]):.2e}"
        print(line)
    return 0


def cmd_evolve(args: argparse.Namespace) -> int:
    rows = conservation_comparison(args.dim, args.t_final, args.steps, args.n_points)
    wanted = [r for r in rows if args.method in ("all", r["method"])]
    print(f"harmonic coherent state, {args.dim}D, t_final={args.t_final}, steps={args.steps}")
    for r in wanted:
        print(
            f"  {r['method']:15s} norm drift {sci(float(r['max_norm_drift']))}  "
            f"rel. energy drift {sci(float(r['relative_energy_drift']))}  "
            f"L2 error {sci(float(r['l2_error_vs_analytic']))}  "
            f"time {num(float(r['seconds']), 3)} s"
        )
    return 0


def cmd_tunneling(args: argparse.Namespace) -> int:
    barrier = RectangularBarrier(args.height, args.width)
    m = int(args.points) | 1
    h = args.width / m
    res = transmission_fd(args.energy, args.height, m, h)
    ex = float(barrier.transmission(args.energy)[0])
    print(f"barrier V0={args.height}, a={m * h:.4f} ({m} grid points), E={args.energy}")
    print(f"  T numerical (FD)  = {res.transmission[0]:.8f}")
    print(f"  T exact           = {ex:.8f}   abs error {abs(res.transmission[0] - ex):.2e}")
    print(f"  R + T - 1         = {res.reflection[0] + res.transmission[0] - 1:.2e}")
    return 0


def cmd_poisson(args: argparse.Namespace) -> int:
    f, u = manufactured_problem(args.n)
    r = solve_poisson(f, args.solver, tol=args.tol)
    print(
        f"{args.solver}, n={args.n}: {r.iterations} iterations, setup {r.setup_time:.3f} s, "
        f"solve {r.solve_time:.3f} s, final relative residual {r.residuals[-1]:.2e}, "
        f"L2 error vs exact {l2_error(r.x, u):.3e}"
    )
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    readme = Path(args.readme) if args.readme else None
    generate_report(Path(args.out), readme, quick=args.quick)
    return 0


def cmd_list(_: argparse.Namespace) -> int:
    print("stationary solvers :", ", ".join(SOLVERS))
    print("time-dependent     :", ", ".join(TIME_METHODS))
    print("poisson solvers    :", ", ".join(POISSON_SOLVERS))
    print("potentials         :", ", ".join(POTENTIALS))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pdeqs", description="Educational PDE / Schrodinger solvers and their convergence."
    )
    p.add_argument("--version", action="version", version=f"pdeqs {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    c = sub.add_parser("convergence", help="empirical order of convergence of one method")
    c.add_argument("method", choices=CONVERGENCE_TARGETS)
    c.add_argument("--potential", choices=[*POTENTIALS, "all"], default="harmonic")
    c.add_argument("--dim", type=int, choices=(1, 2), default=1, help="time-dependent methods only")
    c.add_argument("--quick", action="store_true", help="fewer refinements (fast)")
    c.add_argument("--no-plot", action="store_true")
    c.add_argument("--out", default="figures", help="directory for the PNG")
    c.set_defaults(func=cmd_convergence)

    e = sub.add_parser("eigen", help="lowest eigenvalues of a 1D potential")
    e.add_argument("--potential", choices=list(POTENTIALS), default="harmonic")
    e.add_argument("--method", choices=list(SOLVERS), default="numerov")
    e.add_argument("-n", "--n-points", type=int, default=400)
    e.add_argument("--states", type=int, default=6)
    e.set_defaults(func=cmd_eigen)

    t = sub.add_parser("evolve", help="norm/energy/accuracy of the propagators (coherent state)")
    t.add_argument("--method", choices=["all", *TIME_METHODS], default="all")
    t.add_argument("--dim", type=int, choices=(1, 2), default=1)
    t.add_argument("--t-final", type=float, default=10.0)
    t.add_argument("--steps", type=int, default=5000)
    t.add_argument("-n", "--n-points", type=int, default=None)
    t.set_defaults(func=cmd_evolve)

    u = sub.add_parser("tunneling", help="transmission through a rectangular barrier")
    u.add_argument("--energy", type=float, default=0.72)
    u.add_argument("--height", type=float, default=1.0)
    u.add_argument("--width", type=float, default=1.55)
    u.add_argument("--points", type=int, default=81, help="grid points across the barrier")
    u.set_defaults(func=cmd_tunneling)

    q = sub.add_parser("poisson", help="solve a 2D Poisson problem")
    q.add_argument("--solver", choices=POISSON_SOLVERS, default="mg-cg")
    q.add_argument("-n", type=int, default=127, help="interior points per side (2^k - 1)")
    q.add_argument("--tol", type=float, default=1e-10)
    q.set_defaults(func=cmd_poisson)

    r = sub.add_parser("report", help="regenerate all figures, tables and the README tables")
    r.add_argument("--out", default="docs")
    r.add_argument("--readme", default=None)
    r.add_argument("--quick", action="store_true")
    r.set_defaults(func=cmd_report)

    sub.add_parser("list", help="list available solvers").set_defaults(func=cmd_list)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
