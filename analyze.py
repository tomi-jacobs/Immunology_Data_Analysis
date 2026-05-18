"""
analyze.py
Tomi Jacobs
May 2026

Purpose:
    Analyze immune cell data from teiko.db

    Part 2: compute_cell_frequencies()
        - Calculate percentage of each cell type per sample

    Part 3: compare_responders()
        - Compare cell frequencies between responders and non-responders
        - Melanoma patients on miraclib, PBMC samples only
        - Mann-Whitney U test for statistical significance

    Part 4: analyze_baseline_subset()
        - Explore melanoma miraclib PBMC samples at baseline (time = 0)
"""

import sqlite3
import pandas as pd
from scipy import stats


# -- Database file --
DB_FILE = "teiko.db"


# -- Connect to database --
def get_connection():
    """Return a connection to the SQLite database"""
    conn = sqlite3.connect(DB_FILE)
    return conn


# ── Part 2 ─────────────────────────────────────────────────────────────────────
def compute_cell_frequencies(conn):
    """
    For each sample, calculate what percentage of total cells
    each immune cell population represents.

    Steps:
        1. Get all cell counts from the database
        2. Calculate total cells per sample (sum of all 5 populations)
        3. Convert wide format (5 columns) to long format (1 row per population)
        4. Calculate percentage = count / total * 100

    Returns a DataFrame with columns:
        sample, total_count, population, count, percentage
    """

    # Step 1: Pull cell counts from database
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

    # Step 2: Calculate total cell count per sample
    # sum across all 5 cell type columns for each row
    cell_columns = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    data["total_count"] = data[cell_columns].sum(axis=1)

    # Step 3: Melt from wide to long format
    # Before melt: 1 row per sample, 5 cell type columns
    # After melt:  5 rows per sample, 1 column for population name, 1 for count
    long_data = data.melt(
        id_vars=["sample", "total_count"],
        value_vars=cell_columns,
        var_name="population",
        value_name="count"
    )

    # Step 4: Calculate percentage
    long_data["percentage"] = (long_data["count"] / long_data["total_count"] * 100).round(2)

    # Reorder columns to match required output
    long_data = long_data[["sample", "total_count", "population", "count", "percentage"]]

    # Sort by sample then population
    long_data = long_data.sort_values(["sample", "population"]).reset_index(drop=True)

    return long_data


# ── Part 3 ─────────────────────────────────────────────────────────────────────
def compare_responders(conn):
    """
    Compare cell population frequencies between responders and non-responders.

    Filter: melanoma + miraclib + PBMC + known response (yes or no) only

    Statistical test: Mann-Whitney U test
        - Non-parametric test (does not assume normal distribution)
        - Compares two independent groups (responders vs non-responders)
        - p < 0.05 means the difference is statistically significant

    Returns:
        filtered_data  -- the filtered DataFrame
        stats_results  -- DataFrame with test results per cell population
    """

    # Pull full dataset with all metadata
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

    # Calculate total count and percentage for each cell type
    cell_columns = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
    data["total_count"] = data[cell_columns].sum(axis=1)
    for col in cell_columns:
        data[col + "_pct"] = (data[col] / data["total_count"] * 100).round(2)

    # Filter to melanoma + miraclib + PBMC + known response only
    filtered = data[
        (data["condition"] == "melanoma") &
        (data["treatment"] == "miraclib") &
        (data["sample_type"] == "PBMC") &
        (data["response"].isin(["yes", "no"]))
    ].copy()

    print(f"Part 3 - Filtered samples: {len(filtered)}")
    print(f"  Responders:     {(filtered['response'] == 'yes').sum()}")
    print(f"  Non-responders: {(filtered['response'] == 'no').sum()}")

    # Run Mann-Whitney U test for each cell population
    results = []

    for col in cell_columns:

        pct_col = col + "_pct"

        # Split into two groups
        responders = filtered[filtered["response"] == "yes"][pct_col]
        non_responders = filtered[filtered["response"] == "no"][pct_col]

        # Run the test
        stat, pval = stats.mannwhitneyu(responders, non_responders, alternative="two-sided")

        # Store results
        results.append({
            "population":              col,
            "responder_median_%":      round(responders.median(), 2),
            "non_responder_median_%":  round(non_responders.median(), 2),
            "p_value":                 round(pval, 4),
            "significant (p<0.05)":    "YES" if pval < 0.05 else "no"
        })

    stats_results = pd.DataFrame(results)
    return filtered, stats_results


# ── Part 4 ─────────────────────────────────────────────────────────────────────
def analyze_baseline_subset(conn):
    """
    Explore melanoma PBMC samples at baseline (time = 0) on miraclib.

    Answers 4 questions:
        Q1: How many samples from each project?
        Q2: How many subjects were responders vs non-responders?
        Q3: How many subjects were male vs female?
        Q4: Average B cells for melanoma males, responders, time = 0
    """

    # Filter directly in SQL query for efficiency
    query = """
        SELECT s.sample,
               s.subject,
               s.response,
               s.treatment,
               s.sample_type,
               s.time_from_treatment_start,
               subj.condition,
               subj.sex,
               subj.project,
               cc.b_cell
        FROM samples s
        JOIN cell_counts cc ON s.sample = cc.sample
        JOIN subjects subj ON s.subject = subj.subject
        WHERE subj.condition            = 'melanoma'
          AND s.treatment               = 'miraclib'
          AND s.sample_type             = 'PBMC'
          AND s.time_from_treatment_start = 0
    """
    data = pd.read_sql_query(query, conn)

    print(f"\nPart 4 - Baseline subset: {len(data)} samples")

    # Q1: How many samples from each project?
    q1 = data.groupby("project")["sample"].count().reset_index()
    q1.columns = ["project", "sample_count"]

    # Q2: How many unique subjects were responders vs non-responders?
    # Use nunique() to count patients not samples
    q2 = data.groupby("response")["subject"].nunique().reset_index()
    q2.columns = ["response", "subject_count"]

    # Q3: How many unique subjects were male vs female?
    q3 = data.groupby("sex")["subject"].nunique().reset_index()
    q3.columns = ["sex", "subject_count"]

    # Q4: Average B cells for melanoma males who responded at time = 0
    male_responders = data[(data["sex"] == "M") & (data["response"] == "yes")]
    avg_b_cell = round(male_responders["b_cell"].mean(), 2)

    return {
        "q1_samples_per_project":       q1,
        "q2_subjects_per_response":     q2,
        "q3_subjects_per_sex":          q3,
        "q4_avg_bcell_male_responders": avg_b_cell
    }


# -- Main: run all parts and print results --
def main():

    conn = get_connection()

    # Part 2
    print("=" * 50)
    print("PART 2: Cell Population Frequencies")
    print("=" * 50)
    freq = compute_cell_frequencies(conn)
    print(f"Total rows: {len(freq):,} ({freq['sample'].nunique():,} samples x 5 populations)")
    print()
    print(freq.head(10).to_string(index=False))

    # Part 3
    print()
    print("=" * 50)
    print("PART 3: Responder vs Non-Responder Comparison")
    print("=" * 50)
    filtered, stats_results = compare_responders(conn)
    print()
    print("Mann-Whitney U Test Results:")
    print(stats_results.to_string(index=False))

    # Part 4
    print()
    print("=" * 50)
    print("PART 4: Baseline Subset Analysis")
    print("=" * 50)
    results = analyze_baseline_subset(conn)

    print()
    print("Q1 - Samples per project:")
    print(results["q1_samples_per_project"].to_string(index=False))

    print()
    print("Q2 - Subjects per response:")
    print(results["q2_subjects_per_response"].to_string(index=False))

    print()
    print("Q3 - Subjects per sex:")
    print(results["q3_subjects_per_sex"].to_string(index=False))

    print()
    print(f"Q4 - Average B cells (melanoma males, responders, time=0): "
          f"{results['q4_avg_bcell_male_responders']:.2f}")

    conn.close()


# Run main when script is executed directly
if __name__ == "__main__":
    main()
