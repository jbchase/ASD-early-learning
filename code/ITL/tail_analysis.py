import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from trimming import prepare_itl_dataframe, add_percentile_trim_flags

def _standardize_text(s):
    if pd.isna(s):
        return s
    s = str(s).strip()
    mapping = {
        "shank3b": "Shank3B",
        "tsc2": "TSC2",
        "deterministic": "Deterministic",
        "probabilistic": "Probabilistic",
        "male": "male",
        "female": "female",
        "het": "Het",
        "wt": "WT",
    }
    return mapping.get(s.lower(), s)


def run_tail_analysis(
    input_csv: str,
    output_folder: str,
    percentile: int = 95,
    min_itl: float = 0,
) -> None:
    """
    Reads trial_level_itl_all_rows.csv, computes genotype-blind animal×session trimming,
    and saves a retained-vs-excluded histogram figure on a log10 ITL scale.
    """
    input_path = Path(input_csv)
    output_dir = Path(output_folder)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)

    required = ["A_ID", "session", "genotype", "strain", "schedule_label", "sex", "itl"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input file is missing required columns: {missing}")

    for col in ["genotype", "strain", "schedule_label", "sex"]:
        df[col] = df[col].map(_standardize_text)

    df = prepare_itl_dataframe(
        df,
        min_itl=min_itl,
        drop_na_itl=True,
    )

    df["cohort"] = (
        df["strain"].astype(str)
        + " | "
        + df["sex"].astype(str)
        + " | "
        + df["schedule_label"].astype(str)
    )

    df, animal_session_thresholds = add_percentile_trim_flags(
        df,
        group_cols=["A_ID", "session"],
        percentiles=[percentile],
        value_col="itl",
        prefix="animal_session",
    )

    keep_col = f"keep_animal_session_p{percentile}"
    excl_col = f"excluded_animal_session_p{percentile}"

    # Save enriched trial-level file with the new flags
    df.to_csv(output_dir / f"trial_level_with_animal_session_p{percentile}_flags.csv", index=False)

    # Log-transform for plotting
    df["log10_itl"] = np.log10(df["itl"])

    cohort_order = sorted(df["cohort"].dropna().unique())
    session_order = sorted(df["session"].dropna().unique())

    nrows = len(cohort_order)
    ncols = len(session_order)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5 * ncols, 3.2 * nrows),
        sharex=True,
        sharey=False,
    )

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    elif ncols == 1:
        axes = axes[:, np.newaxis]

    for i, cohort in enumerate(cohort_order):
        for j, session in enumerate(session_order):
            ax = axes[i, j]
            sub = df[(df["cohort"] == cohort) & (df["session"] == session)].copy()

            if sub.empty:
                ax.set_visible(False)
                continue

            retained = sub.loc[sub[keep_col], "log10_itl"].values
            excluded = sub.loc[sub[excl_col], "log10_itl"].values

            bins = np.linspace(sub["log10_itl"].min(), sub["log10_itl"].max(), 40)

            ax.hist(retained, bins=bins, density=True, alpha=0.7, label="Retained")
            ax.hist(excluded, bins=bins, density=True, alpha=0.7, label="Excluded")

            ax.set_title(f"{cohort}\nSession {int(session)}")
            ax.set_xlabel("log10(ITL seconds)", fontsize=9)
            ax.set_ylabel("Density", fontsize=9)
            ax.tick_params(axis="both", labelsize=8)
            ax.set_ylim(0, 1.5)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7)

    fig.suptitle(
        f"Retained vs excluded ITL distributions\n"
        f"(animal×session {percentile}th percentile trim; collapsed across genotype)",
        y=1.01
    )
    fig.tight_layout()

    fig.savefig(
        output_dir / f"retained_vs_excluded_log10_itl_histograms_p{percentile}.png",
        dpi=300,
        bbox_inches="tight",
    )

    fig.savefig(
        output_dir / f"retained_vs_excluded_log10_itl_histograms_p{percentile}.svg",
        bbox_inches="tight",
    )
    plt.close(fig)

    print(f"Tail analysis figure saved to: {output_dir}")
