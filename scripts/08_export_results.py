"""Final export checkpoint. Extend once model figures and feature importance are generated."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
needed = [
    ROOT / "outputs/tables/model_comparison.csv",
    ROOT / "data/reports/antibiotic_eligibility.csv",
]
missing = [p for p in needed if not p.exists()]
if missing:
    raise SystemExit("Missing required outputs:\n" + "\n".join(str(p.relative_to(ROOT)) for p in missing))
print("Core result tables are present. Add final figures/reports before marking M10 complete.")
