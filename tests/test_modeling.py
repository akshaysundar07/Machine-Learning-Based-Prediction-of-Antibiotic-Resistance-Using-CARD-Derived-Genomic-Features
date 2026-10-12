"""Tests for M6-M8: fold discipline, leakage, metric integrity, determinism, output completeness."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kpamr.modeling.cv import (MODEL_NAMES, cross_validate_model, imbalance_ratio, lineage_folds,  # noqa: E402
                               load_antibiotic, load_settings, make_model, metric_row)

SUMMARY = pd.read_csv(ROOT / "data/reports/ml_dataset_summary.csv")
CODES = [f.replace("_ml.csv", "") for f in SUMMARY.filename]
METRICS = ROOT / "models/metrics"


def test_sixteen_eligible_antibiotics():
    assert len(CODES) == 16 and "CST" not in CODES and "TGC" not in CODES


@pytest.mark.parametrize("code", ["GEN", "IPM"])
def test_features_have_no_label_or_metadata_columns(code):
    X, y, ids, fold = load_antibiotic(ROOT, code)
    assert X.shape[1] == 145 and all(c.startswith(("gene__", "variant__")) for c in X.columns)
    assert set(y.unique()) == {0, 1} and sorted(fold.unique()) == [0, 1, 2, 3, 4]


def test_class_weight_rule_uses_training_labels_only():
    balanced = np.array([0] * 50 + [1] * 50); skewed = np.array([0] * 80 + [1] * 20)
    assert imbalance_ratio(skewed) >= load_settings()["imbalance_ratio"] > imbalance_ratio(balanced)
    assert make_model("logistic_regression", skewed).class_weight == "balanced"
    assert make_model("logistic_regression", balanced).class_weight is None
    assert make_model("xgboost", skewed).scale_pos_weight == pytest.approx(4.0)
    assert make_model("xgboost", balanced).scale_pos_weight == 1.0


def test_validation_fold_is_never_in_training(monkeypatch):
    import kpamr.modeling.cv as cv

    X, y, ids, fold = load_antibiotic(ROOT, "SXT")
    seen = []
    fit_batches = []

    model_class = type(make_model("logistic_regression", y))
    original_fit = model_class.fit
    original_choose = cv.choose_threshold

    def fit_spy(self, Xtr, ytr, *args, **kwargs):
        seen.append(set(Xtr.index))
        return original_fit(self, Xtr, ytr, *args, **kwargs)

    def threshold_spy(name, X_train, y_train, seed=42):
        start = len(seen)
        threshold = original_choose(
            name, X_train, y_train, seed
        )
        fit_batches.append((start, len(seen)))
        return threshold

    monkeypatch.setattr(model_class, "fit", fit_spy)
    monkeypatch.setattr(cv, "choose_threshold", threshold_spy)

    cv.cross_validate_model("logistic_regression", X, y, fold)

    assert len(fit_batches) == fold.nunique()

    for k, (start, end) in zip(
        sorted(fold.unique()), fit_batches
    ):
        expected_train = set(X.index[fold != k])
        validation = set(X.index[fold == k])

        # Three inner fits for threshold selection.
        inner_fits = seen[start:end]
        assert len(inner_fits) == 3

        for train_indices in inner_fits:
            assert train_indices < expected_train
            assert train_indices.isdisjoint(validation)

        # The next fit trains the final outer-fold model.
        outer_fit = seen[end]
        assert outer_fit == expected_train
        assert outer_fit.isdisjoint(validation)

    assert len(seen) == fold.nunique() * 4


def test_oof_predictions_cover_every_isolate_once_and_are_deterministic():
    X, y, ids, fold = load_antibiotic(ROOT, "CAZ")
    pf1, o1 = cross_validate_model("random_forest", X, y, fold, 42)
    pf2, o2 = cross_validate_model("random_forest", X, y, fold, 42)
    assert not np.isnan(o1).any() and len(o1) == len(y)
    np.testing.assert_array_equal(o1, o2)
    assert pf1.n_valid.sum() == len(y) and len(pf1) == 5


def test_metric_row_confusion_matrix_and_auc():
    r = metric_row(np.array([0, 0, 1, 1]), np.array([0, 1, 1, 1]), np.array([.1, .6, .7, .9]))
    assert (r["tn"], r["fp"], r["fn"], r["tp"]) == (1, 1, 0, 2)
    assert r["accuracy"] == .75 and r["recall"] == 1.0 and r["roc_auc"] == 1.0


@pytest.mark.parametrize("col", ["MLST", "K locus"])
def test_lineage_folds_never_split_a_group(col):
    X, y, ids, fold = load_antibiotic(ROOT, "GEN")
    f = lineage_folds(ROOT, ids, y, col)
    g = pd.read_csv(ROOT / "data/raw/metadata/Metadata.csv", dtype={"id": str}).set_index("id").loc[ids, col].astype(str).to_numpy()
    assert (f >= 0).all() and all(f[g == v].nunique() == 1 for v in set(g))


@pytest.mark.parametrize("model", [*MODEL_NAMES, "majority_baseline"])
def test_saved_metrics_complete_and_valid(model):
    d = pd.read_csv(METRICS / f"{model}.csv")
    per_fold = d[d.fold != "mean_std"]
    assert per_fold.antibiotic.nunique() == 16 and len(per_fold) == 80
    expected = dict(zip(CODES, SUMMARY.total_samples))
    assert per_fold.groupby("antibiotic").n_valid.sum().to_dict() == expected
    assert per_fold[["accuracy", "precision", "recall", "f1"]].apply(lambda s: s.between(0, 1)).all().all()
    oof = pd.read_csv(METRICS / f"{model}_oof_predictions.csv")
    assert oof.y_prob.between(0, 1).all() and oof.y_prob.notna().all()


def test_baseline_is_chance_and_models_beat_it():
    comp = pd.read_csv(ROOT / "outputs/tables/model_comparison.csv")
    base = comp[comp.model == "majority_baseline"]
    assert (base.roc_auc_mean.dropna() == 0.5).all()
    for m in MODEL_NAMES:
        assert comp[comp.model == m].roc_auc_mean.mean() > 0.8


def test_model_comparison_table_integrity():
    comp = pd.read_csv(ROOT / "outputs/tables/model_comparison.csv")
    assert len(comp) == 16 * 4 and (comp.tn + comp.fp + comp.fn + comp.tp).eq(comp.n_isolates).all()
    assert (comp.n_resistant + comp.n_susceptible).eq(comp.n_isolates).all()
    assert comp[comp.best_model_for_antibiotic].groupby("antibiotic").size().eq(1).all()
    assert comp[comp.best_model_by_f1].groupby("antibiotic").size().eq(1).all()
    best = comp[comp.best_model_for_antibiotic]
    assert best.model.isin(MODEL_NAMES).all()


def test_lineage_robustness_has_all_schemes():
    r = pd.read_csv(ROOT / "outputs/tables/lineage_robustness.csv")
    assert set(r.scheme) == {"standard_stratified", "grouped_MLST", "grouped_K_locus"} and len(r) == 16 * 3 * 3


def test_figures_exist():
    for f in ("roc_curves_all.png", "confusion_matrices_best_model.png", "auc_comparison.png"):
        assert (ROOT / "outputs/figures" / f).stat().st_size > 10_000


def test_coefficients_direction_and_checkpoints():
    c = pd.read_csv(METRICS / "logistic_regression_coefficients.csv")
    assert c.antibiotic.nunique() == 16 and (c.groupby("antibiotic").size() == 145).all()
    assert ((c.coefficient > 0) == (c.direction == "toward_resistant")).all()
    assert len(list((ROOT / "models/checkpoints").glob("*.joblib"))) == 48


def test_settings_come_from_config_files_not_code():
    cfg = load_settings()
    project = __import__("yaml").safe_load((ROOT / "configs/project.yaml").read_text())
    assert cfg["seed"] == project["project"]["random_seed"] and cfg["n_splits"] == project["modeling"]["validation"]["n_splits"]
    # hyperparameters actually reach the estimators
    y = np.array([0] * 50 + [1] * 50)
    assert make_model("random_forest", y).n_estimators == cfg["hyperparameters"]["random_forest"]["n_estimators"]
    assert make_model("xgboost", y).max_depth == cfg["hyperparameters"]["xgboost"]["max_depth"]
    assert make_model("logistic_regression", y).C == cfg["hyperparameters"]["logistic_regression"]["C"]
    assert make_model("logistic_regression", y).random_state == cfg["seed"]

@pytest.mark.parametrize("model", MODEL_NAMES)
def test_saved_threshold_predictions_match_metrics(model):
    oof = pd.read_csv(METRICS / f"{model}_oof_predictions.csv")
    metrics = pd.read_csv(METRICS / f"{model}.csv")

    assert not oof.duplicated(["antibiotic", "isolate_id"]).any()
    assert oof.threshold.between(0.05, 0.95).all()
    assert oof.y_pred.eq(
        (oof.y_prob >= oof.threshold).astype(int)
    ).all()

    for (antibiotic, fold), group in oof.groupby(
        ["antibiotic", "fold"]
    ):
        assert group.threshold.nunique() == 1

        rows = metrics[
            metrics.antibiotic.eq(antibiotic)
            & metrics.fold.astype(str).eq(str(fold))
        ]
        assert len(rows) == 1
        row = rows.iloc[0]

        calculated = metric_row(
            group.y_true,
            group.y_pred,
            group.y_prob,
        )

        assert np.isclose(row.threshold, group.threshold.iloc[0])
        assert row.n_valid == len(group)

        for metric, value in calculated.items():
            assert np.isclose(row[metric], value)


def test_threshold_selection_only_receives_outer_training_data(monkeypatch):
    import kpamr.modeling.cv as cv

    X, y, ids, fold = load_antibiotic(ROOT, "IPM")
    seen = []
    original = cv.choose_threshold

    def spy(name, X_train, y_train, seed=42):
        seen.append(set(X_train.index))
        return original(name, X_train, y_train, seed)

    monkeypatch.setattr(cv, "choose_threshold", spy)
    cv.cross_validate_model("logistic_regression", X, y, fold)

    assert len(seen) == fold.nunique()

    for k, received_indices in zip(sorted(fold.unique()), seen):
        expected = set(X.index[fold != k])
        assert received_indices == expected