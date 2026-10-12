"""M8 evaluation: fair per-antibiotic model comparison, confusion matrices, ROC curves, lineage robustness.

Reads models/metrics/*.csv (from M6/M7). Outputs to outputs/tables and outputs/figures.
Model choice is by ROC-AUC (primary) then F1 - never accuracy alone.
"""
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kpamr.modeling.cv import MODEL_NAMES, cross_validate_model, lineage_folds, load_antibiotic, load_settings, summarize  # noqa: E402

CFG = load_settings()
THRESHOLD, SEED, N_SPLITS = CFG["threshold"], CFG["seed"], CFG["n_splits"]
GRID_COLS = 4


def grid(n, cell=(4, 3.8)):
    """Subplot grid sized from the number of antibiotics; unused panels are hidden."""
    import math
    rows = math.ceil(n / GRID_COLS)
    fig, axes = plt.subplots(rows, GRID_COLS, figsize=(cell[0] * GRID_COLS, cell[1] * rows), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    return fig, axes.ravel()[:n]


MODELS = list(MODEL_NAMES)
mdir = ROOT / "models/metrics"
tables = ROOT / "outputs/tables"; figs = ROOT / "outputs/figures"
tables.mkdir(parents=True, exist_ok=True); figs.mkdir(parents=True, exist_ok=True)
summary = pd.read_csv(ROOT / "data/reports/ml_dataset_summary.csv")
codes = [f.replace("_ml.csv", "") for f in summary.filename]
label = dict(zip(codes, summary.antibiotic))

missing = [m for m in [*MODELS, "majority_baseline"] if not (mdir / f"{m}.csv").exists()]
if missing:
    raise SystemExit(f"Run scripts/06_train_models.py first; missing metrics for {missing}")
folds = pd.concat([pd.read_csv(mdir / f"{m}.csv") for m in [*MODELS, "majority_baseline"]], ignore_index=True)
oof = pd.concat([pd.read_csv(mdir / f"{m}_oof_predictions.csv", dtype={"isolate_id": str}) for m in [*MODELS, "majority_baseline"]], ignore_index=True)

# ---- 1. comparison table + pooled out-of-fold confusion matrices
rows, cms = [], []
for code in codes:
    for m in [*MODELS, "majority_baseline"]:
        g = oof[(oof.antibiotic == code) & (oof.model == m)]
        pred = g.y_pred.astype(int)
        tn, fp, fn, tp = confusion_matrix(g.y_true, pred, labels=[0, 1]).ravel()
        s = folds[(folds.antibiotic == code) & (folds.model == m) & (folds.fold == "mean_std")].iloc[0]
        rows.append({"antibiotic": label[code], "model": m, "n_isolates": len(g), "n_resistant": int(g.y_true.sum()),
                     "n_susceptible": int((g.y_true == 0).sum()), "resistant_fraction": round(g.y_true.mean(), 3),
                     **{c: s[c] for c in s.index if c.endswith(("_mean", "_std"))},
                     "pooled_roc_auc": roc_auc_score(g.y_true, g.y_prob),
                     "tn": tn, "fp": fp, "fn": fn, "tp": tp})
        cms.append({"antibiotic": label[code], "model": m, "TN": tn, "FP": fp, "FN": fn, "TP": tp})
comp = pd.DataFrame(rows)
real = comp[comp.model != "majority_baseline"].copy()
real["rank"] = real.groupby("antibiotic")[["roc_auc_mean", "f1_mean"]].apply(
    lambda d: d.apply(tuple, axis=1).rank(ascending=False, method="first")).reset_index(level=0, drop=True)
comp["best_model_for_antibiotic"] = comp.index.map(lambda i: bool(real["rank"].get(i, 0) == 1))
# second view: best at the configured decision threshold (F1 then ROC-AUC). AUC measures ranking only; F1/recall show
# whether resistant isolates are actually caught. They disagree for the carbapenems, so both are reported.
real["rank_f1"] = real.groupby("antibiotic")[["f1_mean", "roc_auc_mean"]].apply(
    lambda d: d.apply(tuple, axis=1).rank(ascending=False, method="first")).reset_index(level=0, drop=True)
comp["best_model_by_f1"] = comp.index.map(lambda i: bool(real["rank_f1"].get(i, 0) == 1))
comp.round(4).sort_values(["antibiotic", "roc_auc_mean"], ascending=[True, False]).to_csv(tables / "model_comparison.csv", index=False)
pd.DataFrame(cms).to_csv(tables / "confusion_matrices.csv", index=False)

overall = (comp.groupby("model")[[c for c in comp if c.endswith("_mean")]].mean().round(4)
           .join(comp[comp.best_model_for_antibiotic].groupby("model").size().rename("n_antibiotics_best")).fillna(0))
overall.to_csv(tables / "model_overall_summary.csv")

# ---- 2. ROC curves (pooled out-of-fold) per antibiotic + combined grid
colors = {"logistic_regression": "tab:blue", "random_forest": "tab:green", "xgboost": "tab:red"}
fig, axes = grid(len(codes))
for ax, code in zip(axes, codes):
    for m in MODELS:
        g = oof[(oof.antibiotic == code) & (oof.model == m)]
        fpr, tpr, _ = roc_curve(g.y_true, g.y_prob)
        ax.plot(fpr, tpr, color=colors[m], label=f"{m.replace('_', ' ')} ({roc_auc_score(g.y_true, g.y_prob):.2f})")
    ax.plot([0, 1], [0, 1], "k:", lw=1); ax.set_title(f"{label[code]} (n={len(g)})"); ax.legend(fontsize=7, loc="lower right")
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
fig.suptitle(f"Pooled out-of-fold ROC curves ({N_SPLITS}-fold stratified CV, seed {SEED})"); fig.tight_layout()
fig.savefig(figs / "roc_curves_all.png", dpi=130); plt.close(fig)

# ---- 3. confusion matrices of the best model per antibiotic
best = comp[comp.best_model_for_antibiotic].set_index("antibiotic")
fig, axes = grid(len(codes), cell=(3.5, 3.5))
for ax, code in zip(axes, codes):
    b = best.loc[label[code]]; cm = np.array([[b.tn, b.fp], [b.fn, b.tp]])
    ax.imshow(cm, cmap="Blues")
    for (i, j), v in np.ndenumerate(cm): ax.text(j, i, int(v), ha="center", va="center", color="white" if v > cm.max() / 2 else "black")
    ax.set_xticks([0, 1], ["S", "R"]); ax.set_yticks([0, 1], ["S", "R"]); ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
    best_name = b["model"].replace("_", " ")
    ax.set_title(f"{label[code]} - {best_name}", fontsize=9)
fig.suptitle("Confusion matrices (pooled out-of-fold, best model by ROC-AUC then F1)"); fig.tight_layout()
fig.savefig(figs / "confusion_matrices_best_model.png", dpi=130); plt.close(fig)

# ---- 4. AUC comparison bar chart with CV std
fig, ax = plt.subplots(figsize=(max(8, 0.95 * len(codes)), 5)); w = 0.27
for i, m in enumerate(MODELS):
    d = comp[comp.model == m].set_index("antibiotic").loc[[label[c] for c in codes]]
    ax.bar(np.arange(len(codes)) + (i - 1) * w, d.roc_auc_mean, w, yerr=d.roc_auc_std, label=m.replace("_", " "), color=colors[m], capsize=2)
ax.axhline(0.5, color="k", ls=":", lw=1, label="chance (AUC 0.5)"); ax.set_xticks(range(len(codes)), [label[c] for c in codes])
ax.set_ylabel("ROC-AUC (mean ± SD over 5 folds)"); ax.set_ylim(0.3, 1.05); ax.legend(ncol=4); fig.tight_layout()
fig.savefig(figs / "auc_comparison.png", dpi=130); plt.close(fig)

# ---- 5. lineage-aware robustness (grouped CV by MLST and by K locus)
schemes = {"standard_stratified": None, **{f"grouped_{c.replace(' ', '_')}": c for c in CFG["lineage_columns"]}}
rob = []
for code in codes:
    X, y, ids, fold = load_antibiotic(ROOT, code)
    for scheme, group_col in schemes.items():
        f = fold if group_col is None else lineage_folds(ROOT, ids, y, group_col)
        for m in MODELS:
            pf, _ = cross_validate_model(m, X, y, f, SEED)
            rob.append({"antibiotic": label[code], "model": m, "scheme": scheme, "n_folds": f.nunique(), **summarize(pf)})
rob = pd.DataFrame(rob)
rob.round(4).to_csv(tables / "lineage_robustness.csv", index=False)
piv = rob.pivot_table(index="model", columns="scheme", values="roc_auc_mean").round(4)
for sch in schemes:
    if sch != "standard_stratified":
        piv[f"drop_{sch}_vs_standard"] = (piv.standard_stratified - piv[sch]).round(4)
piv.to_csv(tables / "lineage_robustness_summary.csv")
print(overall.to_string()); print(piv.to_string())
print("best by F1 counts:", comp[comp.best_model_by_f1].model.value_counts().to_dict())
print("best model counts:", comp[comp.best_model_for_antibiotic].model.value_counts().to_dict())
