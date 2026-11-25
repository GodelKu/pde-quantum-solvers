"""Repeated wall-clock measurements: one warm-up run, then the median and range of N runs."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TypeVar

import numpy as np

T = TypeVar("T")

#: Timed repetitions behind every wall time in the README (after one untimed warm-up run).
REPEATS = 5


def repeat_runs(run: Callable[[], T], repeats: int = REPEATS, warmup: bool = True) -> list[T]:
    """Call ``run`` once untimed (if ``warmup``), then ``repeats`` times; return those results."""
    if warmup:
        run()
    return [run() for _ in range(repeats)]


def time_stats(seconds: Sequence[float]) -> dict[str, float]:
    """Median, minimum and maximum of repeated wall times (seconds)."""
    s = np.asarray(seconds, dtype=float)
    return {
        "median": float(np.median(s)),
        "min": float(s.min()),
        "max": float(s.max()),
        "repeats": float(len(s)),
    }
