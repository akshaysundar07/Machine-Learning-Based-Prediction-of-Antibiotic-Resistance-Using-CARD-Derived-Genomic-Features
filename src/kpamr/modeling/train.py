from dataclasses import dataclass
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier


@dataclass(frozen=True)
class ModelSpec:
    name: str
    estimator: object


def default_models(seed: int = 42) -> list[ModelSpec]:
    return [
        ModelSpec(
            "logistic_regression",
            LogisticRegression(max_iter=5000, class_weight="balanced", random_state=seed),
        ),
        ModelSpec(
            "random_forest",
            RandomForestClassifier(
                n_estimators=500,
                class_weight="balanced",
                random_state=seed,
                n_jobs=-1,
            ),
        ),
        ModelSpec(
            "xgboost",
            XGBClassifier(
                n_estimators=300,
                learning_rate=0.05,
                max_depth=3,
                subsample=0.8,
                colsample_bytree=0.8,
                eval_metric="logloss",
                random_state=seed,
                n_jobs=-1,
            ),
        ),
    ]
