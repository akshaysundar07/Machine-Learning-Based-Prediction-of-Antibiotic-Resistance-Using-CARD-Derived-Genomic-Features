import pandas as pd


def assert_unique_key(df: pd.DataFrame, key: str, name: str = "table") -> None:
    dup = df[df.duplicated(key, keep=False)]
    if not dup.empty:
        raise ValueError(f"{name} has duplicate {key} values: {dup[key].nunique()} duplicated keys")


def class_counts(df: pd.DataFrame, antibiotic_col="antibiotic", phenotype_col="phenotype_std") -> pd.DataFrame:
    return (
        df.groupby([antibiotic_col, phenotype_col], dropna=False)
        .size()
        .rename("n")
        .reset_index()
        .sort_values([antibiotic_col, phenotype_col])
    )
