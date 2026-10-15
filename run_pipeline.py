"""Run the original RobuProt scripts from the repository root.

WARNING: Research prototype. Outputs are not clinically validated.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
STEPS = ["step1.py", "step2.py", "step3.py", "step4.py", "step3.py", "save_results.py", "report.py"]

def main():
    for folder in ("data", "results", "figures"):
        (ROOT / folder).mkdir(exist_ok=True)
    # Remove previous generated files so step3 does not accidentally use stale scores.
    generated = [
        ROOT / "data" / name for name in ("intensity_matrix.csv", "metadata.csv", "intensity_norm.csv", "intensity_batch_corrected.csv")
    ]
    generated += list((ROOT / "results").glob("*.csv")) + list((ROOT / "results").glob("*.txt"))
    generated += list((ROOT / "figures").glob("*.png"))
    generated += [ROOT / "proteomics_report.docx"]
    for path in generated:
        if path.exists():
            path.unlink()
    for step in STEPS:
        print(f"\n=== Running {step} ===", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts" / step)], cwd=ROOT, check=True)
    print("\nPipeline finished; review scientific limitations in README.md.")

if __name__ == "__main__":
    main()
