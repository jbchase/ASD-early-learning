
import pandas as pd
from typing import Sequence


def prepare_itl_dataframe(
    df: pd.DataFrame,
    min_itl: float | None = None,
    drop_na_itl: bool = True,
) -> pd.DataFrame:
    """
    Standardize key fields used for ITL trimming.

    - coerces session and itl to numeric
    - strips A_ID
    - optionally drops NaN ITLs
    - optionally filters ITLs <= min_itl
    """
    out = df.copy()

    required = ["A_ID", "session", "itl"]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"Input dataframe is missing required columns: {missing}")

    out["A_ID"] = out["A_ID"].astype(str).str.strip()
    out["session"] = pd.to_numeric(out["session"], errors="coerce")
    out["itl"] = pd.to_numeric(out["itl"], errors="coerce")

    if drop_na_itl:
        out = out.loc[out["itl"].notna()].copy()

    if min_itl is not None:
        out = out.loc[out["itl"] > min_itl].copy()

    return out


def add_percentile_trim_flags(
    df: pd.DataFrame,
    group_cols: Sequence[str],
    percentiles: Sequence[int],
    value_col: str = "itl",
    prefix: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Compute percentile thresholds within groups and append keep/excluded flags.

    Returns
    -------
    enriched_df : pd.DataFrame
        Original dataframe with threshold + keep/excluded columns added.
    thresholds : pd.DataFrame
        One row per group with percentile thresholds.
    """
    if not group_cols:
        raise ValueError("group_cols must contain at least one column")
    if value_col not in df.columns:
        raise ValueError(f"'{value_col}' column not found in dataframe")

    missing = [c for c in group_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Grouping columns missing from dataframe: {missing}")

    if prefix is None:
        prefix = "_".join(group_cols)

    q_values = [p / 100.0 for p in percentiles]

    thresholds = (
        df.groupby(list(group_cols))[value_col]
        .quantile(q_values)
        .unstack()
        .rename(columns={p / 100.0: f"{prefix}_threshold_p{p}" for p in percentiles})
        .reset_index()
    )

    enriched = df.merge(thresholds, on=list(group_cols), how="left")

    for p in percentiles:
        threshold_col = f"{prefix}_threshold_p{p}"
        keep_col = f"keep_{prefix}_p{p}"
        excl_col = f"excluded_{prefix}_p{p}"

        enriched[keep_col] = enriched[value_col] <= enriched[threshold_col]
        enriched[excl_col] = ~enriched[keep_col]

    return enriched, thresholds
