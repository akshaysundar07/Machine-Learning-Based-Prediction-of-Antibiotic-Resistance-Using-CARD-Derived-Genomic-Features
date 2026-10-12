from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kpamr.modeling.cv import (
    MODEL_NAMES,
    load_antibiotic,
    make_model,
    load_settings,
)


def main():
    dictionary = pd.read_csv(
        ROOT / "data/reports/amr_feature_dictionary.csv"
    ).set_index("feature_id")

    summary = pd.read_csv(
        ROOT / "data/reports/ml_dataset_summary.csv"
    )
    codes = [
        filename.replace("_ml.csv", "")
        for filename in summary.filename
    ]

    seed = load_settings()["seed"]
    rows = []

    for code in codes:
        X, y, ids, folds = load_antibiotic(ROOT, code)
        features = X.columns.tolist()

        for model_name in MODEL_NAMES:
            for fold in sorted(folds.unique()):
                train = folds.ne(fold)
                X_train = X.loc[train]
                y_train = y.loc[train]

                model = make_model(
                    model_name, y_train, seed
                ).fit(X_train, y_train)

                if model_name == "logistic_regression":
                    values = model.coef_[0]
                    scores = np.abs(values)

                elif model_name == "random_forest":
                    values = model.feature_importances_
                    scores = values

                else:
                    gains = model.get_booster().get_score(
                        importance_type="gain"
                    )
                    values = np.array([
                        gains.get(feature, 0.0)
                        for feature in features
                    ])
                    scores = values

                frame = pd.DataFrame({
                    "antibiotic": code,
                    "model": model_name,
                    "fold": int(fold),
                    "feature_id": features,
                    "display_name": dictionary.loc[
                        features, "display_name"
                    ].to_numpy(),
                    "value": values,
                    "ranking_score": scores,
                    "training_present_count": (
                        X_train.sum().to_numpy()
                    ),
                })

                frame = frame.sort_values(
                    ["ranking_score", "feature_id"],
                    ascending=[False, True],
                )
                frame["rank"] = np.arange(1, len(frame) + 1)

                # Do not count zero-score features as important.
                frame["top10"] = (
                    frame["rank"].le(10)
                    & frame.ranking_score.gt(0)
                )
                frame["positive"] = frame.value.gt(0)
                frame["negative"] = frame.value.lt(0)

                rows.append(frame)

            print(f"{code}: {model_name} finished")

    details = pd.concat(rows, ignore_index=True)

    stability = details.groupby(
        ["antibiotic", "model", "feature_id", "display_name"],
        as_index=False,
    ).agg(
        top10_folds=("top10", "sum"),
        mean_rank=("rank", "mean"),
        rank_std=("rank", "std"),
        positive_folds=("positive", "sum"),
        negative_folds=("negative", "sum"),
        minimum_training_presence=("training_present_count", "min"),
    )

    # Direction counts apply only to LR coefficients.
    tree_rows = stability.model.ne("logistic_regression")
    stability.loc[
        tree_rows, ["positive_folds", "negative_folds"]
    ] = np.nan

    stability = stability.sort_values(
        ["antibiotic", "model", "top10_folds",
         "mean_rank", "feature_id"],
        ascending=[True, True, False, True, True],
    )

    destination = ROOT / "outputs/tables"
    destination.mkdir(parents=True, exist_ok=True)

    details.to_csv(
        destination / "feature_importance_by_fold.csv",
        index=False,
    )
    stability.to_csv(
        destination / "feature_stability.csv",
        index=False,
    )

    print("\nIPM features appearing in the top 10 in at least 4/5 folds:")
    print(
        stability[
            stability.antibiotic.eq("IPM")
            & stability.top10_folds.ge(4)
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()