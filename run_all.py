"""Run every analysis script in order, then re-test every claim.

    python run_all.py        (about 2 minutes on a laptop)
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for script in sorted((ROOT / "analysis").glob("0*.py")):
    print(f"\n=== {script.name} ===", flush=True)
    subprocess.run([sys.executable, str(script)], check=True, cwd=ROOT)
print("\n=== check_claims.py ===", flush=True)
sys.exit(subprocess.run([sys.executable, str(ROOT / "check_claims.py")], cwd=ROOT).returncode)
