"""M6-M8 core: cross-validation on the saved M5 folds.

Rules enforced here:
- fold_id == k is validation; all other folds are training (never the reverse).
- Nothing is tuned: hyperparameters are fixed, so there is no tuning on validation data.
- Class weighting is decided from the TRAINING fold only, and only when imbalance warrants it.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import yaml
import pandas as pd
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, average_precision_score, confusion_matrix,
                             f1_score, precision_score, recall_score, roc_auc_score)
from xgboost import XGBClassifier

MODEL_NAMES = ("logistic_regression", "random_forest", "xgboost")
META_COLS = ("isolate_id", "antibiotic", "binary_label")
ROOT = Path(__file__).resolve().parents[3]


@lru_cache(maxsize=None)
def load_settings(root: Path = ROOT) -> dict:
    """Read seed / n_splits from project.yaml and everything else from modeling.yaml. Nothing is hard-coded."""
    project = yaml.safe_load((root / "configs/project.yaml").read_text(encoding="utf-8"))
    extra = yaml.safe_load((root / "configs/modeling.yaml").read_text(encoding="utf-8"))
    m = extra["modeling"]
    return {"seed": int(project["project"]["random_seed"]),
            "n_splits": int(project["modeling"]["validation"]["n_splits"]),
            "threshold": float(m["decision_threshold"]),
            "imbalance_ratio": float(m["class_weight_min_imbalance_ratio"]),
            "lineage_columns": list(m["lineage_group_columns"]),
            "hyperparameters": m["hyperparameters"],
            "expected": extra["audit"]["expected"]}


def imbalance_ratio(y) -> float:
    pos = int(np.sum(y == 1)); neg = int(np.sum(y == 0))
    return max(pos, neg) / max(min(pos, neg), 1)


def make_model(name: str, y_train, seed: int | None = None):
    cfg = load_settings()
    seed = cfg["seed"] if seed is None else seed
    weighted = imbalance_ratio(y_train) >= cfg["imbalance_ratio"]
    cw = "balanced" if weighted else None
    if name == "majority_baseline":
        return DummyClassifier(strategy="prior")
    hp = dict(cfg["hyperparameters"][name])
    if name == "logistic_regression":
        return LogisticRegression(**hp, class_weight=cw, random_state=seed)
    if name == "random_forest":
        return RandomForestClassifier(**hp, class_weight=cw, random_state=seed, n_jobs=1)
    if name == "xgboost":
        pos = int(np.sum(y_train == 1)); neg = int(np.sum(y_train == 0))
        return XGBClassifier(**hp, eval_metric="logloss", random_state=seed, n_jobs=1,
                             scale_pos_weight=(neg / max(pos, 1)) if weighted else 1.0)
    raise ValueError(name)


def load_antibiotic(root: Path, code: str):
    """Return (X, y, ids, fold_id) after verifying dataset and folds agree exactly."""
    d = pd.read_csv(root / f"data/processed/ml_datasets/{code}_ml.csv", dtype={"isolate_id": str})
    f = pd.read_csv(root / f"data/processed/splits/{code}_folds.csv", dtype={"isolate_id": str})
    if not (d.isolate_id.tolist() == f.isolate_id.tolist() and d.binary_label.tolist() == f.binary_label.tolist()):
        raise ValueError(f"{code}: ML dataset and fold file disagree")
    X = d.drop(columns=list(META_COLS))
    return X, d.binary_label.astype(int), d.isolate_id, f.fold_id.astype(int)


def metric_row(y_true, y_pred, y_prob) -> dict:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    both = len(np.unique(y_true)) == 2
    return {"accuracy": accuracy_score(y_true, y_pred),
            "precision": precision_score(y_true, y_pred, zero_division=0),
            "recall": recall_score(y_true, y_pred, zero_division=0),
            "f1": f1_score(y_true, y_pred, zero_division=0),
            "roc_auc": roc_auc_score(y_true, y_prob) if both else np.nan,
            "pr_auc": average_precision_score(y_true, y_prob) if both else np.nan,
            "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)}


def cross_validate_model(name, X, y, fold_id, seed=None):
    """Choose thresholds inside training data; evaluate on saved outer folds."""
    cfg = load_settings()
    seed = cfg["seed"] if seed is None else seed

    rows = []
    oof = np.full(len(y), np.nan)

    for k in sorted(fold_id.unique()):
        va = (fold_id == k).to_numpy()
        tr = ~va

        X_train = X.loc[tr]
        y_train = y.loc[tr]

        threshold = (
            cfg["threshold"]
            if name == "majority_baseline"
            else choose_threshold(name, X_train, y_train, seed)
        )

        model = make_model(name, y_train, seed)
        model.fit(X_train, y_train)

        prob = model.predict_proba(X.loc[va])[:, 1]
        pred = (prob >= threshold).astype(int)
        oof[va] = prob

        rows.append({
            "fold": int(k),
            "threshold": threshold,
            "n_train": int(tr.sum()),
            "n_valid": int(va.sum()),
            **metric_row(y.loc[va].to_numpy(), pred, prob),
        })

    return pd.DataFrame(rows), oof


def summarize(per_fold: pd.DataFrame) -> dict:
    out = {}
    for m in ("accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc"):
        out[f"{m}_mean"] = per_fold[m].mean(); out[f"{m}_std"] = per_fold[m].std(ddof=1)
    return out


def lineage_folds(root: Path, ids: pd.Series, y: pd.Series, group_col: str, seed: int | None = None, n_splits: int | None = None) -> pd.Series:
    """Grouped stratified folds: isolates sharing MLST (or K locus) never straddle train/validation."""
    from sklearn.model_selection import StratifiedGroupKFold
    cfg = load_settings()
    seed = cfg["seed"] if seed is None else seed
    n_splits = cfg["n_splits"] if n_splits is None else n_splits
    meta = pd.read_csv(root / "data/raw/metadata/Metadata.csv", dtype={"id": str}).set_index("id")
    groups = meta.loc[ids, group_col].astype(str).to_numpy()
    fold = np.full(len(y), -1)
    for k, (_, va) in enumerate(StratifiedGroupKFold(n_splits, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups)):
        fold[va] = k
    return pd.Series(fold, index=y.index)

def choose_threshold(name, X_train, y_train, seed=42):
    """Choose a threshold using inner CV on the outer training set."""
    from sklearn.model_selection import StratifiedKFold

    inner_cv = StratifiedKFold(
        n_splits=3,
        shuffle=True,
        random_state=seed,
    )

    inner_prob = np.full(len(y_train), np.nan)

    for train_idx, valid_idx in inner_cv.split(X_train, y_train):
        model = make_model(
            name,
            y_train.iloc[train_idx],
            seed,
        )

        model.fit(
            X_train.iloc[train_idx],
            y_train.iloc[train_idx],
        )

        inner_prob[valid_idx] = model.predict_proba(
            X_train.iloc[valid_idx]
        )[:, 1]

    thresholds = np.linspace(0.05, 0.95, 91)

    scores = [
        f1_score(
            y_train,
            inner_prob >= threshold,
            zero_division=0,
        )
        for threshold in thresholds
    ]

    # When scores tie, prefer the threshold closest to 0.5.
    best_score = max(scores)
    candidates = [
        threshold
        for threshold, score in zip(thresholds, scores)
        if np.isclose(score, best_score, rtol=0, atol=1e-12)
    ]

    return float(min(candidates, key=lambda t: abs(t - 0.5)))
