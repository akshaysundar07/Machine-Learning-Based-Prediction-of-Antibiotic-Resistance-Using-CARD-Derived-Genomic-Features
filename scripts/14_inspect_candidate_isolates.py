from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ["517281", "681054"]

features = pd.read_csv(
    ROOT / "data/processed/master/amr_feature_matrix.csv",
    dtype={"isolate_id": str},
).set_index("isolate_id")

dictionary = pd.read_csv(
    ROOT / "data/reports/amr_feature_dictionary.csv"
).set_index("feature_id")

lookup = pd.read_csv(
    ROOT / "outputs/tables/isolate_literature_lookup.csv",
    dtype={"isolate_id": str},
).set_index("isolate_id")

rows = []

for isolate_id in CANDIDATES:
    info = lookup.loc[isolate_id]

    print(f"\n=== MRSN{isolate_id} ===")
    print(
        info[[
            "biosample", "genbank_accessions",
            "MLST", "IPM", "MEM", "ETP",
        ]].to_string()
    )

    present = features.columns[
        features.loc[isolate_id].eq(1)
    ]

    annotations = dictionary.loc[
        present, ["display_name", "feature_type"]
    ]

    print("\nRecorded genomic determinants:")
    print(annotations.to_string(index=False))

    for feature_id, annotation in annotations.iterrows():
        rows.append({
            "isolate_id": isolate_id,
            "feature_id": feature_id,
            "display_name": annotation.display_name,
            "feature_type": annotation.feature_type,
        })

    print("\nOut-of-fold IPM predictions:")
    for model in [
        "logistic_regression", "random_forest", "xgboost"
    ]:
        predictions = pd.read_csv(
            ROOT / f"models/metrics/{model}_oof_predictions.csv",
            dtype={"isolate_id": str},
        )

        match = predictions[
            predictions.isolate_id.eq(isolate_id)
            & predictions.antibiotic.eq("IPM")
        ]

        if len(match) != 1:
            raise ValueError(
                f"Expected one IPM prediction: {isolate_id}/{model}"
            )

        print(match[[
            "model", "y_true", "y_prob", "threshold", "y_pred"
        ]].to_string(index=False))

destination = ROOT / "outputs/tables/candidate_isolate_features.csv"
pd.DataFrame(rows).to_csv(destination, index=False)
print(f"\nSaved: {destination.name}")