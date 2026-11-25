import json

import numpy as np
import pytest

from schrodinger_bench.report import _clean, _fit_frequency
from schrodinger_bench.tables import sci, timing


def test_frequency_fit_ignores_fast_oscillations():
    """Slow cosine plus a fast one of a quarter of its amplitude, as in the double-well run."""
    t = np.arange(0.0, 700.0, 0.02)
    slow = 0.0097
    y = 0.5 + 0.45 * np.cos(slow * t + 0.1) + 0.1 * np.cos(2.1 * t) + 0.02 * np.cos(0.33 * t)
    w, _, _ = _fit_frequency(t, y, guess=0.0098, filter_width=15.0)
    assert w == pytest.approx(slow, rel=1e-6)


def test_results_json_is_standard_json():
    data = _clean({"a": np.float64(np.nan), "b": [np.inf, -np.inf, np.int64(3)], "c": (1.5,)})
    text = json.dumps(data, allow_nan=False)
    assert json.loads(text) == {"a": None, "b": ["inf", "-inf", 3], "c": [1.5]}
    assert sci(json.loads(text)["b"][0]) == "inf" and sci(None) == "-"


def test_timing_cell_shows_median_and_range():
    assert timing({"median": 0.25, "min": 0.2, "max": 0.5}, digits=2) == "0.25 (0.20-0.50)"
