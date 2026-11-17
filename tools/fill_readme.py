"""Fill the README tables from docs/results/results.json without rerunning the studies.

Every block between ``<!-- BEGIN:name -->`` and ``<!-- END:name -->`` is rendered by
``pdeqs.tables`` from the JSON alone, and ``tests/test_readme_sync.py`` fails if README,
``docs/results/tables.md`` and the JSON disagree.  ``pdeqs report --out docs --readme README.md``
reruns everything; this script only refreshes the README after the prose around the tables was
edited.  Prose outside the blocks avoids result numbers.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pdeqs.tables import fill, load_results, render_tables, tables_markdown  # noqa: E402


def main() -> int:
    tables = render_tables(load_results(ROOT / "docs" / "results" / "results.json"))
    readme = ROOT / "README.md"
    text = readme.read_text(encoding="utf-8")
    new, found = fill(text, tables)
    missing = sorted(set(tables) - set(found))
    if missing:
        print(f"README has no marker for: {', '.join(missing)}")
        return 1
    readme.write_text(new, encoding="utf-8")
    (ROOT / "docs" / "results" / "tables.md").write_text(tables_markdown(tables), encoding="utf-8")
    print("README updated" if new != text else "README already up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
