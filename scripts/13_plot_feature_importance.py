from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = pd.read_csv(
        ROOT / "outputs/tables/feature_importance.csv"
    )

    destination = ROOT / "outputs/figures/feature_importance"
    destination.mkdir(parents=True, exist_ok=True)

    models = [
        "logistic_regression",
        "random_forest",
        "xgboost",
    ]

    axis_labels = {
        "logistic_regression": "Signed coefficient",
        "random_forest": "Mean decrease in impurity",
        "xgboost": "Average split gain",
    }

    for antibiotic in sorted(data.antibiotic.unique()):
        fig, axes = plt.subplots(
            1, 3, figsize=(18, 6), constrained_layout=True
        )

        for ax, model in zip(axes, models):
            top = data[
                data.antibiotic.eq(antibiotic)
                & data.model.eq(model)
                & data["rank"].le(10)
                & data.ranking_score.gt(0)
            ].sort_values("rank", ascending=False)

            colors = (
                [
                    "tab:red" if value > 0 else "tab:blue"
                    for value in top.value
                ]
                if model == "logistic_regression"
                else "tab:green"
            )

            ax.barh(top.display_name, top.value, color=colors)
            ax.axvline(0, color="black", linewidth=0.8)
            ax.set_title(model.replace("_", " ").title())
            ax.set_xlabel(axis_labels[model])

        label = "C/T" if antibiotic == "CT" else antibiotic
        fig.suptitle(
            f"{label}: top genomic predictors\n"
            "Full-data fitted models; score scales differ across models"
        )

        fig.savefig(
            destination / f"{antibiotic}_top_features.png",
            dpi=180,
        )
        plt.close(fig)

    print("Saved 16 plots in outputs/figures/feature_importance")
    print("LR: red = toward resistant; blue = toward susceptible")
    print("Tree importance does not indicate association direction")


if __name__ == "__main__":
    main()