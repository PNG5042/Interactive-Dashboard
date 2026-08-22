import sqlite3

import matplotlib
matplotlib.use("Agg")  # headless-safe: needed for environments with no display (e.g. Codespaces)
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DB_FILE = "cell_counts.db"
BOXPLOT_FILE = "response_boxplots.png"
STATS_FILE = "response_stats.csv"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def load_subset(conn):
    """melanoma + miraclib + PBMC samples, with relative frequency per population."""
    query = """
        SELECT
            sm.sample_id AS sample,
            sm.subject_id,
            s.response,
            cc.population,
            cc.count
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        JOIN cell_counts cc ON cc.sample_id = sm.sample_id
        WHERE s.condition = 'melanoma'
          AND s.treatment = 'miraclib'
          AND sm.sample_type = 'PBMC'
          AND s.response IS NOT NULL
    """
    df = pd.read_sql_query(query, conn)

    totals = df.groupby("sample")["count"].sum().rename("total_count")
    df = df.join(totals, on="sample")
    df["percentage"] = df["count"] / df["total_count"] * 100
    return df


def make_boxplots(df):
    fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharey=False)

    for ax, pop in zip(axes, POPULATIONS):
        subset = df[df["population"] == pop]
        sns.boxplot(data=subset, x="response", y="percentage", order=["no", "yes"], ax=ax)
        sns.stripplot(data=subset, x="response", y="percentage", order=["no", "yes"],
                       ax=ax, color="black", alpha=0.25, size=2, jitter=0.2)
        ax.set_title(pop)
        ax.set_xlabel("responder")
        ax.set_ylabel("% of total cells" if pop == POPULATIONS[0] else "")

    fig.suptitle("Miraclib melanoma PBMC samples: responders (yes) vs non-responders (no)")
    fig.tight_layout()
    fig.savefig(BOXPLOT_FILE, dpi=150)
    print(f"saved boxplots -> {BOXPLOT_FILE}")


def run_stats(df):
    subject_level = (
        df.groupby(["subject_id", "response", "population"])["percentage"]
        .mean()
        .reset_index()
    )

    rows = []
    for pop in POPULATIONS:
        subset = subject_level[subject_level["population"] == pop]
        yes = subset.loc[subset["response"] == "yes", "percentage"]
        no = subset.loc[subset["response"] == "no", "percentage"]

        stat, p = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "n_responder_subjects": len(yes),
            "n_non_responder_subjects": len(no),
            "median_responder": yes.median(),
            "median_non_responder": no.median(),
            "u_stat": stat,
            "p_value": p,
        })

    results = pd.DataFrame(rows)

    # correct for the fact we're running 5 tests
    reject, p_adj, _, _ = multipletests(results["p_value"], method="fdr_bh")
    results["p_adj"] = p_adj
    results["significant"] = reject

    results = results.sort_values("p_value").reset_index(drop=True)
    return results


if __name__ == "__main__":
    conn = sqlite3.connect(DB_FILE)
    df = load_subset(conn)
    conn.close()

    print(f"samples in subset: {df['sample'].nunique()}")
    print(f"responders: {df.loc[df['response'] == 'yes', 'sample'].nunique()}, "
          f"non-responders: {df.loc[df['response'] == 'no', 'sample'].nunique()}")
    print()

    make_boxplots(df)

    results = run_stats(df)
    results.to_csv(STATS_FILE, index=False)

    print()
    print(results.to_string(index=False))
    print()
    sig = results.loc[results["significant"], "population"].tolist()
    if sig:
        print(f"significant populations (BH-adjusted p < 0.05): {sig}")
    else:
        print("no populations reached significance after correction")