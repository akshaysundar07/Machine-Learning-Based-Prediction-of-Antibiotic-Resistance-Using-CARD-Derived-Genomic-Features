import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

metadata = pd.read_csv(
    ROOT / "data/raw/metadata/Metadata.csv",
    dtype={"id": str},
)
metadata["id"] = metadata["id"].str.strip()

assert metadata["id"].is_unique, "Duplicate metadata IDs"

records = []

with (ROOT / "data/assembly_data_report.jsonl").open(
    encoding="utf-8"
) as stream:
    for line in stream:
        if not line.strip():
            continue

        record = json.loads(line)
        info = record.get("assemblyInfo", {})
        biosample = info.get("biosample", {})

        names = [
            biosample.get("strain", ""),
            record.get("organism", {})
            .get("infraspecificNames", {})
            .get("strain", ""),
        ]

        names.extend(
            item.get("value", "")
            for item in biosample.get("sampleIds", [])
            if item.get("label") == "Sample name"
        )

        ids = set()
        for name in names:
            match = re.fullmatch(
                r"MRSN[\s_-]*(\d+)", str(name).strip(), re.I
            )
            if match:
                ids.add(match.group(1))

        if len(ids) != 1:
            raise ValueError(f"Ambiguous isolate identity: {names}")

        accession = record.get("accession", "")
        paired = (
            record.get("pairedAccession")
            or info.get("pairedAssembly", {}).get("accession", "")
        )

        records.append({
            "isolate_id": ids.pop(),
            "biosample": biosample.get("accession", ""),
            "accessions": [accession, paired],
        })

rows = []

for isolate_id, group in pd.DataFrame(records).groupby("isolate_id"):
    biosamples = sorted(set(group.biosample) - {""})
    if len(biosamples) != 1:
        raise ValueError(f"Conflicting BioSamples for {isolate_id}")

    accessions = sorted({
        accession
        for values in group.accessions
        for accession in values
        if accession
    })

    rows.append({
        "isolate_id": isolate_id,
        "strain": f"MRSN{isolate_id}",
        "biosample": biosamples[0],
        "genbank_accessions": ";".join(
            a for a in accessions if a.startswith("GCA_")
        ),
        "refseq_accessions": ";".join(
            a for a in accessions if a.startswith("GCF_")
        ),
    })

mapping = pd.DataFrame(rows)

lookup = metadata.merge(
    mapping,
    left_on="id",
    right_on="isolate_id",
    how="left",
    validate="one_to_one",
)

if lookup.biosample.isna().any():
    missing = lookup.loc[lookup.biosample.isna(), "id"].tolist()
    raise ValueError(f"Missing accession mappings: {missing}")

columns = [
    "isolate_id", "strain", "biosample",
    "genbank_accessions", "refseq_accessions",
    "MLST", "Region", "Sample Type",
]

antibiotics = [
    "AMK", "GEN", "TOB", "SAM", "ATM", "CAZ", "CRO",
    "FEP", "CZA", "C/T", "ETP", "IPM", "MEM",
    "CIP", "LVX", "TZP", "SXT", "TET", "TGC",
]

lookup = lookup[columns + antibiotics]

destination = ROOT / "outputs/tables/isolate_literature_lookup.csv"
destination.parent.mkdir(parents=True, exist_ok=True)
lookup.to_csv(destination, index=False)

print(f"Saved lookup for {len(lookup)} isolates")
print(lookup[columns].head().to_string(index=False))