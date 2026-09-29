"""Create human-readable model comparison tables from saved CV metrics."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
src = ROOT / "models/metrics/cv_model_comparison.csv"
if not src.exists():
    raise SystemExit("Run scripts/06_train_models.py first.")
df = pd.read_csv(src)
out = ROOT / "outputs/tables/model_comparison.csv"
out.parent.mkdir(parents=True, exist_ok=True)
df.sort_values(["antibiotic", "roc_auc_mean"], ascending=[True, False]).to_csv(out, index=False)
print(f"Wrote {out.relative_to(ROOT)}")
