"""M6 (Logistic Regression) and M7 (Random Forest, XGBoost) on the M5 datasets + saved folds.

  python scripts/06_train_models.py --stage m6     # logistic regression only
  python scripts/06_train_models.py --stage m7     # random forest + xgboost
  python scripts/06_train_models.py                # all (m6 then m7)

Outputs: models/metrics/<model>.csv (per-fold + summary), models/metrics/<model>_oof_predictions.csv,
models/metrics/logistic_regression_coefficients.csv, models/checkpoints/<ABX>_<model>.joblib (+ .json config).
"""
from pathlib import Path
import argparse, json, sys
import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kpamr.modeling.cv import cross_validate_model, imbalance_ratio, load_antibiotic, load_settings, make_model, summarize  # noqa: E402

SEED = load_settings()["seed"]  # from configs/project.yaml
ap = argparse.ArgumentParser(); ap.add_argument("--stage", choices=["m6", "m7", "all"], default="all")
stage = ap.parse_args().stage
models = {"m6": ["logistic_regression"], "m7": ["random_forest", "xgboost"],
          "all": ["logistic_regression", "random_forest", "xgboost"]}[stage]

summary = pd.read_csv(ROOT / "data/reports/ml_dataset_summary.csv")
codes = [f.replace("_ml.csv", "") for f in summary.filename]
dictionary = pd.read_csv(ROOT / "data/reports/amr_feature_dictionary.csv").set_index("feature_id")
mdir = ROOT / "models/metrics"; cdir = ROOT / "models/checkpoints"
mdir.mkdir(parents=True, exist_ok=True); cdir.mkdir(parents=True, exist_ok=True)

for name in [*models, *(["majority_baseline"] if stage in ("m6", "all") else [])]:
    rows, oofs, coefs = [], [], []
    for code in codes:
        X, y, ids, fold = load_antibiotic(ROOT, code)
        pf, oof = cross_validate_model(name, X, y, fold, SEED)
        pf.insert(0, "antibiotic", code); pf.insert(1, "model", name)
        rows.append(pf)
        thresholds = fold.map(pf.set_index("fold")["threshold"])

        oofs.append(pd.DataFrame({
            "antibiotic": code,
            "model": name,
            "isolate_id": ids,
            "fold": fold,
            "y_true": y,
            "y_prob": oof,
            "threshold": thresholds,
            "y_pred": (oof >= thresholds.to_numpy()).astype(int),
        }))
        if name == "majority_baseline":
            continue
        final = make_model(name, y, SEED).fit(X, y)  # fit on all usable isolates: for interpretation only
        joblib.dump(final, cdir / f"{code.replace('/', '')}_{name}.joblib")
        (cdir / f"{code.replace('/', '')}_{name}.json").write_text(json.dumps(
            {"antibiotic": code, "model": name, "seed": SEED, "n": int(len(y)), "n_features": int(X.shape[1]),
             "imbalance_ratio": round(imbalance_ratio(y), 3), "params": {k: str(v) for k, v in final.get_params().items()}},
            indent=2))
        if name == "logistic_regression":
            c = pd.DataFrame({"antibiotic": code, "feature_id": X.columns, "coefficient": final.coef_[0]})
            c["display_name"] = c.feature_id.map(dictionary.display_name)
            c["direction"] = np.where(c.coefficient > 0, "toward_resistant", np.where(c.coefficient < 0, "toward_susceptible", "none"))
            coefs.append(c)
        print(f"{name:20s} {code:4s} done")
    pf_all = pd.concat(rows, ignore_index=True)
    summ = []
    for code, g in pf_all.groupby("antibiotic", sort=False):
        summ.append({"antibiotic": code, "model": name, "fold": "mean_std", "n_valid": int(g.n_valid.sum()), **summarize(g)})
    pd.concat([pf_all, pd.DataFrame(summ)], ignore_index=True).to_csv(mdir / f"{name}.csv", index=False)
    pd.concat(oofs, ignore_index=True).to_csv(mdir / f"{name}_oof_predictions.csv", index=False)
    if coefs:
        pd.concat(coefs, ignore_index=True).to_csv(mdir / "logistic_regression_coefficients.csv", index=False)
    print("wrote", f"models/metrics/{name}.csv")
