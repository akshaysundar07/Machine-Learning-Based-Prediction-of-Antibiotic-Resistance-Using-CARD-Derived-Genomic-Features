"""Cross-validate Logistic Regression, Random Forest and XGBoost for every eligible antibiotic."""
from pathlib import Path
import sys
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_validate

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kpamr.modeling.train import default_models  # noqa: E402

files = sorted((ROOT / "data/processed/per_antibiotic").glob("*.csv"))
if not files:
    raise SystemExit("No per-antibiotic datasets found. Complete M5 first.")

scoring = {"accuracy": "accuracy", "precision": "precision", "recall": "recall", "f1": "f1", "roc_auc": "roc_auc"}
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
rows = []
for path in files:
    df = pd.read_csv(path)
    y = df.pop("label_binary").astype(int)
    ids = df.pop("isolate_id")
    X = df
    for spec in default_models(seed=42):
        scores = cross_validate(spec.estimator, X, y, cv=cv, scoring=scoring, n_jobs=1)
        row = {"antibiotic": path.stem, "model": spec.name, "n": len(y), "n_features": X.shape[1]}
        for metric in scoring:
            vals = scores[f"test_{metric}"]
            row[f"{metric}_mean"] = vals.mean()
            row[f"{metric}_std"] = vals.std(ddof=1)
        rows.append(row)

out = ROOT / "models/metrics/cv_model_comparison.csv"
out.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out, index=False)
print(f"Wrote {out.relative_to(ROOT)}")
