from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

CASES = {
    "517281": ["blaNDM-5", "blaOXA-181"],
    "681054": ["blaNDM-7"],
}
ANTIBIOTICS = ["IPM", "MEM", "ETP"]
MODELS = ["logistic_regression", "random_forest", "xgboost"]

dictionary = pd.read_csv(
    ROOT / "data/reports/amr_feature_dictionary.csv"
)

carrier_rows = []
prediction_rows = []

for antibiotic in ANTIBIOTICS:
    data = pd.read_csv(
        ROOT / f"data/processed/ml_datasets/{antibiotic}_ml.csv",
        dtype={"isolate_id": str},
    )
    folds = pd.read_csv(
        ROOT / f"data/processed/splits/{antibiotic}_folds.csv",
        dtype={"isolate_id": str},
    )

    data = data.merge(
        folds[["isolate_id", "fold_id"]],
        on="isolate_id",
        validate="one_to_one",
    )

    for isolate_id, genes in CASES.items():
        case = data[data.isolate_id.eq(isolate_id)]
        if len(case) != 1:
            raise ValueError(
                f"Expected one row: {antibiotic}/{isolate_id}"
            )

        case = case.iloc[0]
        training = data[data.fold_id.ne(case.fold_id)]

        for gene in genes:
            matches = dictionary.loc[
                dictionary.display_name.eq(gene), "feature_id"
            ]
            if len(matches) != 1:
                raise ValueError(f"Ambiguous feature: {gene}")

            feature = matches.iloc[0]
            if case[feature] != 1:
                raise ValueError(
                    f"{isolate_id} does not carry {gene}"
                )

            carriers = training[training[feature].eq(1)]

            carrier_rows.append({
                "antibiotic": antibiotic,
                "isolate_id": isolate_id,
                "gene": gene,
                "fold": int(case.fold_id),
                "total_carriers": int(data[feature].sum()),
                "training_carriers": len(carriers),
                "training_resistant": int(
                    carriers.binary_label.eq(1).sum()
                ),
                "training_susceptible": int(
                    carriers.binary_label.eq(0).sum()
                ),
            })

        for model in MODELS:
            oof = pd.read_csv(
                ROOT / f"models/metrics/{model}_oof_predictions.csv",
                dtype={"isolate_id": str},
            )
            match = oof[
                oof.antibiotic.eq(antibiotic)
                & oof.isolate_id.eq(isolate_id)
            ]

            if len(match) != 1:
                raise ValueError(
                    f"Missing/duplicate prediction: "
                    f"{antibiotic}/{isolate_id}/{model}"
                )

            result = match.iloc[0]

            if int(result.y_true) != int(case.binary_label):
                raise ValueError("Prediction target disagrees with dataset")
            if int(result.fold) != int(case.fold_id):
                raise ValueError("Prediction fold disagrees with dataset")

            fixed_prediction = int(result.y_prob >= 0.5)
            selected_prediction = int(result.y_pred)

            prediction_rows.append({
                "antibiotic": antibiotic,
                "isolate_id": isolate_id,
                "model": model,
                "actual_label": int(result.y_true),
                "probability": result.y_prob,
                "threshold": result.threshold,
                "fixed_prediction": fixed_prediction,
                "selected_prediction": selected_prediction,
                "fixed_correct": (
                    fixed_prediction == int(result.y_true)
                ),
                "selected_correct": (
                    selected_prediction == int(result.y_true)
                ),
            })

carriers = pd.DataFrame(carrier_rows)
predictions = pd.DataFrame(prediction_rows)

destination = ROOT / "outputs/tables"
destination.mkdir(parents=True, exist_ok=True)

carriers.to_csv(
    destination / "rare_allele_training_coverage.csv",
    index=False,
)
predictions.to_csv(
    destination / "rare_allele_case_predictions.csv",
    index=False,
)

print("TRAINING COVERAGE")
print(carriers.to_string(index=False))

print("\nCASE PREDICTIONS: 0 = susceptible, 1 = resistant")
print(predictions.round(4).to_string(index=False))