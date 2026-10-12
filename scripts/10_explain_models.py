from pathlib import Path
import sys

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from kpamr.modeling.cv import MODEL_NAMES, load_antibiotic


def main():
    dictionary = pd.read_csv(
        ROOT / "data/reports/amr_feature_dictionary.csv"
    )

    if dictionary.feature_id.duplicated().any():
        raise ValueError("Duplicate feature IDs in dictionary")

    dictionary = dictionary.set_index("feature_id")

    summary = pd.read_csv(
        ROOT / "data/reports/ml_dataset_summary.csv"
    )
    codes = [
        filename.replace("_ml.csv", "")
        for filename in summary.filename
    ]

    results = []

    for code in codes:
        X, y, ids, folds = load_antibiotic(ROOT, code)
        features = X.columns.tolist()

        missing = set(features) - set(dictionary.index)
        if missing:
            raise ValueError(f"Missing biological mappings: {missing}")

        for model_name in MODEL_NAMES:
            path = ROOT / (
                f"models/checkpoints/{code}_{model_name}.joblib"
            )
            if not path.exists():
                raise FileNotFoundError(path)

            model = joblib.load(path)

            if model.n_features_in_ != len(features):
                raise ValueError(f"{code}/{model_name}: feature count mismatch")

            if hasattr(model, "feature_names_in_"):
                if list(model.feature_names_in_) != features:
                    raise ValueError(
                        f"{code}/{model_name}: feature order mismatch"
                    )

            if model_name == "logistic_regression":
                values = model.coef_[0]
                ranking_scores = np.abs(values)
                method = "signed_coefficient"
                directions = np.where(
                    values > 0,
                    "toward_resistant",
                    np.where(values < 0, "toward_susceptible", "none"),
                )

            elif model_name == "random_forest":
                values = model.feature_importances_
                ranking_scores = values
                method = "mean_decrease_in_impurity"
                directions = ["not_available"] * len(features)

            else:
                booster = model.get_booster()

                if booster.feature_names != features:
                    raise ValueError(
                        f"{code}/xgboost: feature names mismatch"
                    )

                gains = booster.get_score(importance_type="gain")
                values = np.array([
                    gains.get(feature, 0.0) for feature in features
                ])
                ranking_scores = values
                method = "average_split_gain"
                directions = ["not_available"] * len(features)

            if not np.isfinite(values).all():
                raise ValueError(f"{code}/{model_name}: invalid scores")

            frame = pd.DataFrame({
                "antibiotic": code,
                "model": model_name,
                "feature_id": features,
                "display_name": dictionary.loc[
                    features, "display_name"
                ].to_numpy(),
                "feature_type": dictionary.loc[
                    features, "feature_type"
                ].to_numpy(),
                "importance_method": method,
                "value": values,
                "ranking_score": ranking_scores,
                "direction": directions,
                "present_count": X.sum().to_numpy(),
                "n_isolates": len(X),
                "analysis_scope": "full_dataset_fitted_model",
            })

            # Feature ID breaks score ties reproducibly.
            frame = frame.sort_values(
                ["ranking_score", "feature_id"],
                ascending=[False, True],
            ).reset_index(drop=True)

            frame["rank"] = np.arange(1, len(frame) + 1)
            results.append(frame)

    combined = pd.concat(results, ignore_index=True)
    top_features = combined[
        combined["rank"].le(10)
        & combined["ranking_score"].gt(0)
    ].copy()

    destination = ROOT / "outputs/tables"
    destination.mkdir(parents=True, exist_ok=True)

    combined.to_csv(
        destination / "feature_importance.csv",
        index=False,
    )
    top_features.to_csv(
        destination / "top10_features.csv",
        index=False,
    )

    print(f"Antibiotics: {len(codes)}")
    print(f"Model-feature rows: {len(combined)}")
    print("Saved feature_importance.csv and top10_features.csv")

    print("\nIPM: top 5 features per model")
    print(
        combined[
            combined.antibiotic.eq("IPM")
            & combined["rank"].le(5)
        ][[
            "model", "rank", "display_name",
            "value", "direction", "present_count",
        ]].to_string(index=False)
    )


if __name__ == "__main__":
    main()