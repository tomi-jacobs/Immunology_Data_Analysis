"""
dashboard.py
Tomi Jacobs
May 2026

Interactive dashboard for Bob Loblaw's clinical trial immune cell analysis.
Displays results from Parts 2, 3, and 4 of the analysis.

Run with:
    streamlit run dashboard.py
"""
import os
import subprocess

# If database doesn't exist, run load_data.py to create it
if not os.path.exists("teiko.db"):
    subprocess.run(["python", "load_data.py"], check=True)
import streamlit as st
import sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# -- Page configuration --
st.set_page_config(
    page_title="Teiko Immune Cell Dashboard",
    page_icon="🧬",
    layout="wide"
)

# -- Database connection --
DB_FILE = "teiko.db"

def get_connection():
    """Return a connection to the SQLite database"""
    conn = sqlite3.connect(DB_FILE)
    return conn


# ── Part 2: Cell Frequencies ───────────────────────────────────────────────────
def compute_cell_frequencies(conn):
    """Calculate percentage of each cell type per sample"""

    query = """
        SELECT s.sample,
               cc.b_cell,
               cc.cd8_t_cell,
               cc.cd4_t_cell,
               cc.nk_cell,
               cc.monocyte
        FROM samples s
        JOIN cell_counts cc ON s.sample = cc.sample
        ORDER BY s.sample
    """
    data = pd.read_sql_query(query, conn)

    cell_columns = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    data["total_count"] = data[cell_columns].sum(axis=1)

    long_data = data.melt(
        id_vars=["sample", "total_count"],
        value_vars=cell_columns,
        var_name="population",
        value_name="count"
    )

    long_data["percentage"] = (long_data["count"] / long_data["total_count"] * 100).round(2)
    long_data = long_data[["sample", "total_count", "population", "count", "percentage"]]
    long_data = long_data.sort_values(["sample", "population"]).reset_index(drop=True)

    return long_data


# ── Part 3: Responder Comparison ───────────────────────────────────────────────
def compare_responders(conn):
    """Compare cell frequencies between responders and non-responders"""

    query = """
        SELECT s.sample,
               s.subject,
               s.response,
               s.treatment,
               s.sample_type,
               subj.condition,
               cc.b_cell,
               cc.cd8_t_cell,
               cc.cd4_t_cell,
               cc.nk_cell,
               cc.monocyte
        FROM samples s
        JOIN cell_counts cc ON s.sample = cc.sample
        JOIN subjects subj ON s.subject = subj.subject
    """
    data = pd.read_sql_query(query, conn)

    cell_columns = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    data["total_count"] = data[cell_columns].sum(axis=1)
    for col in cell_columns:
        data[col + "_pct"] = (data[col] / data["total_count"] * 100).round(2)

    # Filter
    filtered = data[
        (data["condition"] == "melanoma") &
        (data["treatment"] == "miraclib") &
        (data["sample_type"] == "PBMC") &
        (data["response"].isin(["yes", "no"]))
    ].copy()

    # Statistical tests
    results = []
    for col in cell_columns:
        pct_col = col + "_pct"
        responders = filtered[filtered["response"] == "yes"][pct_col]
        non_responders = filtered[filtered["response"] == "no"][pct_col]
        stat, pval = stats.mannwhitneyu(responders, non_responders, alternative="two-sided")
        results.append({
            "population":             col,
            "responder_median_%":     round(responders.median(), 2),
            "non_responder_median_%": round(non_responders.median(), 2),
            "p_value":                round(pval, 4),
            "significant (p<0.05)":   "YES ✓" if pval < 0.05 else "no"
        })

    stats_df = pd.DataFrame(results)
    return filtered, stats_df


# ── Part 4: Baseline Subset ────────────────────────────────────────────────────
def analyze_baseline_subset(conn):
    """Analyze melanoma PBMC baseline samples on miraclib"""

    query = """
        SELECT s.sample,
               s.subject,
               s.response,
               subj.sex,
               subj.project,
               cc.b_cell
        FROM samples s
        JOIN cell_counts cc ON s.sample = cc.sample
        JOIN subjects subj ON s.subject = subj.subject
        WHERE subj.condition              = 'melanoma'
          AND s.treatment                 = 'miraclib'
          AND s.sample_type               = 'PBMC'
          AND s.time_from_treatment_start = 0
    """
    data = pd.read_sql_query(query, conn)

    q1 = data.groupby("project")["sample"].count().reset_index()
    q1.columns = ["project", "sample_count"]

    q2 = data.groupby("response")["subject"].nunique().reset_index()
    q2.columns = ["response", "subject_count"]

    q3 = data.groupby("sex")["subject"].nunique().reset_index()
    q3.columns = ["sex", "subject_count"]

    male_responders = data[(data["sex"] == "M") & (data["response"] == "yes")]
    avg_b_cell = round(male_responders["b_cell"].mean(), 2)

    return {"q1": q1, "q2": q2, "q3": q3, "q4": avg_b_cell}


# ── Dashboard Layout ───────────────────────────────────────────────────────────
def main():

    # Header
    st.markdown("""
        <div style='background:#19325A; padding:24px 32px; border-radius:10px; margin-bottom:24px'>
            <h1 style='color:white; margin:0; font-size:28px'>🧬 Teiko Immune Cell Dashboard</h1>
            <p style='color:#BBBBEE; margin:6px 0 0; font-size:15px'>
                Clinical Trial Immune Monitoring · Bob Loblaw, Loblaw Bio
            </p>
        </div>
    """, unsafe_allow_html=True)

    # Connect to database
    conn = get_connection()

    # Load all data
    freq_df = compute_cell_frequencies(conn)
    filtered_df, stats_df = compare_responders(conn)
    baseline = analyze_baseline_subset(conn)

    conn.close()

    # ── Tabs ──────────────────────────────────────────────────────────────────
    tab1, tab2, tab3 = st.tabs([
        "📊 Part 2 — Cell Frequencies",
        "🔬 Part 3 — Responder Analysis",
        "📋 Part 4 — Baseline Subset"
    ])

    # ── TAB 1: Cell Frequencies ───────────────────────────────────────────────
    with tab1:
        st.subheader("Cell Population Frequencies Per Sample")
        st.markdown("""
            For each sample, the table below shows the relative frequency (percentage)
            of each immune cell population as a proportion of total cells counted.
        """)

        # Summary stats
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Rows", f"{len(freq_df):,}")
        col2.metric("Unique Samples", f"{freq_df['sample'].nunique():,}")
        col3.metric("Cell Populations", "5")

        st.markdown("---")

        # Sample selector
        sample_options = sorted(freq_df["sample"].unique())
        selected_sample = st.selectbox(
            "Select a sample to inspect:",
            sample_options,
            index=0
        )

        sample_data = freq_df[freq_df["sample"] == selected_sample]

        col_a, col_b = st.columns([1, 1])

        with col_a:
            st.markdown(f"**Frequencies for {selected_sample}**")
            st.dataframe(
                sample_data[["population", "count", "percentage"]],
                use_container_width=True,
                hide_index=True
            )

        with col_b:
            # Pie chart for selected sample
            fig, ax = plt.subplots(figsize=(5, 4))
            colors = ["#19325A", "#2E75B6", "#00695C", "#E65100", "#F9A825"]
            ax.pie(
                sample_data["percentage"],
                labels=sample_data["population"],
                autopct="%1.1f%%",
                colors=colors,
                startangle=90
            )
            ax.set_title(f"Cell distribution — {selected_sample}", fontsize=11)
            st.pyplot(fig)
            plt.close()

        st.markdown("---")
        st.markdown("**Full Frequency Table (first 100 rows)**")
        st.dataframe(freq_df.head(100), use_container_width=True, hide_index=True)

    # ── TAB 2: Responder Analysis ─────────────────────────────────────────────
    with tab2:
        st.subheader("Responder vs Non-Responder Comparison")
        st.markdown("""
            Comparing cell population frequencies between **responders** and **non-responders**
            in melanoma patients receiving **miraclib** (PBMC samples only).
            Statistical significance tested using the **Mann-Whitney U test** (p < 0.05).
        """)

        # Filter summary
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Samples", f"{len(filtered_df):,}")
        col2.metric("Responders", f"{(filtered_df['response']=='yes').sum():,}")
        col3.metric("Non-Responders", f"{(filtered_df['response']=='no').sum():,}")

        st.markdown("---")

        # Stats table
        st.markdown("**Statistical Test Results**")

        # Highlight significant rows
        def highlight_significant(row):
            if "YES" in str(row["significant (p<0.05)"]):
                return ["background-color: #E8F5E9"] * len(row)
            return [""] * len(row)

        st.dataframe(
            stats_df.style.apply(highlight_significant, axis=1),
            use_container_width=True,
            hide_index=True
        )

        st.markdown("""
            > **Key finding:** CD4 T cells show a statistically significant difference
            > between responders and non-responders (p=0.0134).
            > Responders have a higher median CD4 T cell frequency (30.22%) compared
            > to non-responders (29.66%), suggesting CD4 T cells may serve as a
            > predictive biomarker for miraclib response.
        """)

        st.markdown("---")

        # Boxplots
        st.markdown("**Boxplots — Responder vs Non-Responder by Cell Population**")

        cell_columns = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

        fig, axes = plt.subplots(1, 5, figsize=(18, 5))
        fig.suptitle(
            "Cell Population Frequencies: Responders vs Non-Responders\n"
            "Melanoma · Miraclib · PBMC",
            fontsize=13, fontweight="bold", y=1.02
        )

        colors = {"yes": "#2E75B6", "no": "#E65100"}
        labels = {"yes": "Responder", "no": "Non-Responder"}

        for i, col in enumerate(cell_columns):
            ax = axes[i]
            pct_col = col + "_pct"

            resp_data = filtered_df[filtered_df["response"] == "yes"][pct_col]
            non_resp_data = filtered_df[filtered_df["response"] == "no"][pct_col]

            bp = ax.boxplot(
                [resp_data, non_resp_data],
                patch_artist=True,
                widths=0.5,
                medianprops=dict(color="white", linewidth=2)
            )

            # Color the boxes
            bp["boxes"][0].set_facecolor("#2E75B6")
            bp["boxes"][1].set_facecolor("#E65100")

            # Get p-value for this population
            p_row = stats_df[stats_df["population"] == col].iloc[0]
            pval = p_row["p_value"]
            sig = "✓ p=" + str(pval) if "YES" in str(p_row["significant (p<0.05)"]) else "p=" + str(pval)

            ax.set_title(col.replace("_", " ").title(), fontsize=10, fontweight="bold")
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Resp", "Non-Resp"], fontsize=9)
            ax.set_ylabel("Frequency (%)", fontsize=9)
            ax.annotate(sig, xy=(0.5, 0.97), xycoords="axes fraction",
                       ha="center", va="top", fontsize=8.5,
                       color="#2E7D32" if "✓" in sig else "#888888")

        plt.tight_layout()
        st.pyplot(fig)
        plt.close()

    # ── TAB 3: Baseline Subset ────────────────────────────────────────────────
    with tab3:
        st.subheader("Baseline Subset Analysis")
        st.markdown("""
            Exploring **melanoma patients on miraclib** using **PBMC samples at baseline (time = 0)**.
        """)

        # Q4 highlight
        st.metric(
            "Q4 — Avg B Cells (Melanoma Males, Responders, Time=0)",
            f"{baseline['q4']:,.2f}"
        )

        st.markdown("---")

        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("**Q1 — Samples per Project**")
            st.dataframe(baseline["q1"], use_container_width=True, hide_index=True)

            # Bar chart
            fig, ax = plt.subplots(figsize=(4, 3))
            ax.bar(baseline["q1"]["project"], baseline["q1"]["sample_count"],
                   color=["#19325A", "#2E75B6"])
            ax.set_title("Samples per Project", fontsize=11)
            ax.set_xlabel("Project")
            ax.set_ylabel("Sample Count")
            for bar in ax.patches:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 2,
                        str(int(bar.get_height())),
                        ha="center", fontsize=10)
            st.pyplot(fig)
            plt.close()

        with col2:
            st.markdown("**Q2 — Subjects per Response**")
            st.dataframe(baseline["q2"], use_container_width=True, hide_index=True)

            fig, ax = plt.subplots(figsize=(4, 3))
            ax.bar(baseline["q2"]["response"], baseline["q2"]["subject_count"],
                   color=["#E65100", "#2E75B6"])
            ax.set_title("Subjects by Response", fontsize=11)
            ax.set_xlabel("Response")
            ax.set_ylabel("Subject Count")
            for bar in ax.patches:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 2,
                        str(int(bar.get_height())),
                        ha="center", fontsize=10)
            st.pyplot(fig)
            plt.close()

        with col3:
            st.markdown("**Q3 — Subjects per Sex**")
            st.dataframe(baseline["q3"], use_container_width=True, hide_index=True)

            fig, ax = plt.subplots(figsize=(4, 3))
            ax.bar(baseline["q3"]["sex"], baseline["q3"]["subject_count"],
                   color=["#C2185B", "#19325A"])
            ax.set_title("Subjects by Sex", fontsize=11)
            ax.set_xlabel("Sex")
            ax.set_ylabel("Subject Count")
            for bar in ax.patches:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 2,
                        str(int(bar.get_height())),
                        ha="center", fontsize=10)
            st.pyplot(fig)
            plt.close()

    # Footer
    st.markdown("---")
    st.markdown("""
        <div style='text-align:center; color:#888; font-size:12px'>
            Tomi Jacobs · PhD Candidate, Computational Biology · University of Illinois Chicago ·
            github.com/tomi-jacobs
        </div>
    """, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
