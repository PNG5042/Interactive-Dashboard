import sqlite3

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import streamlit as st
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

DB_FILE = "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

st.set_page_config(page_title="Loblaw Bio - Cell Population Dashboard", layout="wide")


@st.cache_data
def load_all_data():
    conn = sqlite3.connect(DB_FILE)
    counts = pd.read_sql_query(
        """
        SELECT
            cc.sample_id AS sample,
            cc.population,
            cc.count,
            sm.sample_type,
            sm.time_from_treatment_start,
            s.subject_id,
            s.project_id,
            s.condition,
            s.age,
            s.sex,
            s.treatment,
            s.response
        FROM cell_counts cc
        JOIN samples sm ON sm.sample_id = cc.sample_id
        JOIN subjects s ON s.subject_id = sm.subject_id
        """,
        conn,
    )
    conn.close()

    totals = counts.groupby("sample")["count"].sum().rename("total_count")
    counts = counts.join(totals, on="sample")
    counts["percentage"] = counts["count"] / counts["total_count"] * 100
    return counts


df = load_all_data()

st.title("Loblaw Bio: Immune Cell Population Dashboard")
st.caption("Miraclib clinical trial - immune cell population analysis")

tab1, tab2, tab3 = st.tabs([
    "Part 2: Sample Frequencies",
    "Part 3: Responder vs Non-Responder",
    "Part 4: Baseline Subset",
])

# ---------------------------------------------------------------------------
# Part 2
# ---------------------------------------------------------------------------
with tab1:
    st.header("Relative frequency of each cell population, per sample")

    projects = sorted(df["project_id"].unique())
    selected_project = st.selectbox("Filter by project (optional)", ["All"] + projects)

    sample_view = df[["sample", "total_count", "population", "count", "percentage", "project_id"]].copy()
    if selected_project != "All":
        sample_view = sample_view[sample_view["project_id"] == selected_project]

    sample_view = sample_view.drop(columns="project_id").sort_values(["sample", "population"])
    st.dataframe(sample_view, use_container_width=True, height=400)
    st.download_button(
        "Download full frequency table (CSV)",
        sample_view.to_csv(index=False),
        file_name="cell_frequencies.csv",
    )

# ---------------------------------------------------------------------------
# Part 3
# ---------------------------------------------------------------------------
with tab2:
    st.header("Miraclib melanoma PBMC samples: responders vs non-responders")

    subset = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
        & (df["response"].notna())
    ]

    n_subjects = subset["subject_id"].nunique()
    n_yes = subset.loc[subset["response"] == "yes", "subject_id"].nunique()
    n_no = subset.loc[subset["response"] == "no", "subject_id"].nunique()
    st.write(f"**{n_subjects} subjects** ({n_yes} responders, {n_no} non-responders), "
             f"{subset['sample'].nunique()} PBMC samples total")

    fig, axes = plt.subplots(1, 5, figsize=(20, 5))
    for ax, pop in zip(axes, POPULATIONS):
        pop_data = subset[subset["population"] == pop]
        sns.boxplot(data=pop_data, x="response", y="percentage", order=["no", "yes"], ax=ax)
        ax.set_title(pop)
        ax.set_xlabel("responder")
        ax.set_ylabel("% of total cells" if pop == POPULATIONS[0] else "")
    fig.tight_layout()
    st.pyplot(fig)

    # stats - averaged per subject first to avoid pseudoreplication across timepoints
    st.subheader("Statistical test (Mann-Whitney U, subject-level, BH-corrected)")
    subject_level = (
        subset.groupby(["subject_id", "response", "population"])["percentage"]
        .mean()
        .reset_index()
    )
    rows = []
    for pop in POPULATIONS:
        pop_data = subject_level[subject_level["population"] == pop]
        yes = pop_data.loc[pop_data["response"] == "yes", "percentage"]
        no = pop_data.loc[pop_data["response"] == "no", "percentage"]
        stat, p = mannwhitneyu(yes, no, alternative="two-sided")
        rows.append({
            "population": pop,
            "median_responder": round(yes.median(), 2),
            "median_non_responder": round(no.median(), 2),
            "p_value": p,
        })
    stats_df = pd.DataFrame(rows)
    reject, p_adj, _, _ = multipletests(stats_df["p_value"], method="fdr_bh")
    stats_df["p_adj"] = p_adj
    stats_df["significant"] = reject
    stats_df = stats_df.sort_values("p_value").reset_index(drop=True)

    st.dataframe(stats_df, use_container_width=True)

    if stats_df["significant"].any():
        sig_pops = stats_df.loc[stats_df["significant"], "population"].tolist()
        st.success(f"Significant after correction: {', '.join(sig_pops)}")
    else:
        st.info("No population reaches significance after correcting for multiple testing "
                "(cd4_t_cell trends closest, raw p ≈ 0.01, adjusted p ≈ 0.06).")

# ---------------------------------------------------------------------------
# Part 4
# ---------------------------------------------------------------------------
with tab3:
    st.header("Baseline (t=0) melanoma / miraclib / PBMC subset")

    baseline = df[
        (df["condition"] == "melanoma")
        & (df["treatment"] == "miraclib")
        & (df["sample_type"] == "PBMC")
        & (df["time_from_treatment_start"] == 0)
    ].drop_duplicates("subject_id")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total subjects at baseline", len(baseline))
    with col2:
        st.metric("Responders", int((baseline["response"] == "yes").sum()))
        st.metric("Non-responders", int((baseline["response"] == "no").sum()))
    with col3:
        st.metric("Male", int((baseline["sex"] == "M").sum()))
        st.metric("Female", int((baseline["sex"] == "F").sum()))

    st.subheader("Samples per project")
    st.bar_chart(baseline["project_id"].value_counts())

    st.divider()
    st.subheader("Melanoma male responders at t=0: average B cell count (all sample/treatment types)")
    b_cell_subset = df[
        (df["condition"] == "melanoma")
        & (df["sex"] == "M")
        & (df["response"] == "yes")
        & (df["time_from_treatment_start"] == 0)
        & (df["population"] == "b_cell")
    ]
    avg_b_cell = b_cell_subset["count"].mean()
    st.metric("Average B cell count", f"{avg_b_cell:.2f}", help=f"n = {len(b_cell_subset)} samples")