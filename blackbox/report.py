"""Small logger: every analysis script prints claim-tagged lines to the console AND to
results/<script>.txt, and stores machine-readable numbers in results/metrics/<script>.json.
check_claims.py reads the JSON files and tests every claim in REPORT.md."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
RESULTS.mkdir(exist_ok=True)
(RESULTS / "metrics").mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)


def _clean(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.ndarray):
        return _clean(x.tolist())
    return x


class Log:
    def __init__(self, name):
        self.name = name
        self.lines = []
        self.metrics = {}

    def __call__(self, *parts):
        s = " ".join(str(p) for p in parts)
        print(s)
        sys.stdout.flush()
        self.lines.append(s)

    def claim(self, cid, text, **values):
        """Record a claim line like '[C07] text' and store its numbers under metrics[cid]."""
        self(f"[{cid}] {text}")
        self.metrics.setdefault(cid, {}).update(_clean(values))

    def save(self):
        (RESULTS / f"{self.name}.txt").write_text("\n".join(self.lines) + "\n")
        (RESULTS / "metrics" / f"{self.name}.json").write_text(json.dumps(self.metrics, indent=1))
