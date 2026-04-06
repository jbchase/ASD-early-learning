import pandas as pd
import numpy as np
from pathlib import Path
from trimming import prepare_itl_dataframe, add_percentile_trim_flags

PERCENTILES = [93, 95, 98]

def compute_itl_for_file(csv_path: Path, meta_row: pd.Series) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = pd.read_csv(csv_path)
    df = df.sort_values("trial").reset_index(drop=True)

    trial = pd.to_numeric(df["trial"], errors="coerce")
    center_in = pd.to_numeric(df["center_in"], errors="coerce")
    last_side_out = pd.to_numeric(df["last_side_out"], errors="coerce")

    out = pd.DataFrame({
        "file_name": csv_path.name,
        "A_ID": meta_row["A_ID"],
        "date": meta_row["date"],
        "session": meta_row["session"],
        "genotype": meta_row["genotype"],
        "strain": meta_row["strain"],
        "schedule_label": meta_row["schedule_label"],
        "sex": meta_row["sex"],
        "trial_prev": trial.iloc[:-1].values,
        "trial_next": trial.shift(-1).iloc[:-1].values,
        "last_side_out_prev_trial": last_side_out.iloc[:-1].values,
        "next_center_in": center_in.shift(-1).iloc[:-1].values,
        "itl": (center_in.shift(-1) - last_side_out).iloc[:-1].values,
        "reward_prev_trial": df["reward"].iloc[:-1].values,
    })

    # Context columns from the previous trial
    optional_columns = {
        "outcome": "outcome_prev_trial",
        "trial_types": "trial_type_prev_trial",
        "port_side": "port_side_prev_trial",
        "odor_name": "odor_name_prev_trial",
        "schedule": "schedule_code_prev_trial",
    }
    for source_col, out_col in optional_columns.items():
        if source_col in df.columns:
            out[out_col] = df[source_col].iloc[:-1].values

    return out


def run_analysis(metadata_path: str, csv_folder: str, output_folder: str) -> None:
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(metadata_path, dtype=str)
    meta["file_name"] = meta["file_name"].astype(str)
    meta["A_ID"] = meta["A_ID"].astype(str).str.strip()
    meta["date"] = meta["date"].astype(str)
    meta["session"] = meta["session"].astype(str)

    trial_frames = []
    for csv_path in sorted(Path(csv_folder).glob("*.csv")):
        match = meta.loc[meta["file_name"] == csv_path.name]
        if match.empty:
            raise ValueError(f"No metadata row found for {csv_path.name}")
        trial_frames.append(compute_itl_for_file(csv_path, match.iloc[0]))

    trial_level = pd.concat(trial_frames, ignore_index=True)
    trial_level.to_csv(output_dir / "trial_level_itl_all_rows.csv", index=False)

    valid = prepare_itl_dataframe(
        trial_level,
        min_itl=0,
        drop_na_itl=True,
    )

    valid, session_thresholds = add_percentile_trim_flags(
        valid,
        group_cols=["session"],
        percentiles=PERCENTILES,
        value_col="itl",
        prefix="session",
    )

    valid, animal_session_thresholds = add_percentile_trim_flags(
        valid,
        group_cols=["A_ID", "session"],
        percentiles=PERCENTILES,
        value_col="itl",
        prefix="animal_session",
    )

    session_thresholds["cut_scope"] = "session"
    animal_session_thresholds["cut_scope"] = "animal_session"

    all_thresholds = pd.concat(
        [session_thresholds, animal_session_thresholds],
        ignore_index=True,
        sort=False
    )

    valid.to_csv(output_dir / "trial_level_with_filter_flags.csv", index=False)
    all_thresholds.to_csv(output_dir / "cutoffs_all_groupings.csv", index=False)
  
    # Animal-level summaries for Prism
    summary_rows = []
    for cut_scope, keep_prefix in [
        ("session", "keep_session_p"),
        ("animal_session", "keep_animal_session_p"),
    ]:
        for p in PERCENTILES:
            sub = valid[valid[f"{keep_prefix}{p}"]].copy()

            animal = (
                sub.groupby(["strain", "sex", "schedule_label", "session", "genotype", "A_ID"])
                .agg(
                    n_trials=("itl", "size"),
                    mean_itl=("itl", "mean"),
                    median_itl=("itl", "median"),
                    sem_itl=("itl", lambda s: s.std(ddof=1) / np.sqrt(len(s)) if len(s) > 1 else np.nan),
                    sd_itl=("itl", "std"),
                    min_itl=("itl", "min"),
                    max_itl=("itl", "max"),
                )
                .reset_index()
            )

            animal["cut_scope"] = cut_scope
            animal["percentile"] = p
            summary_rows.append(animal)

    animal_summary_long = pd.concat(summary_rows, ignore_index=True)
    animal_summary_long.to_csv(output_dir / "animal_summary_all_options.csv", index=False)

    # Main analysis now uses animal×session 95th percentile
    main_analysis = animal_summary_long[
        (animal_summary_long["cut_scope"] == "animal_session") &
        (animal_summary_long["percentile"] == 95)
    ].copy()

    main_analysis.to_csv(output_dir / "animal_summary_animal_session_p95.csv", index=False)
