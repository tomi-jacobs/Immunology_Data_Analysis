# Teiko Technical Assessment
**Tomi Jacobs** · PhD Candidate, Computational Biology · University of Illinois Chicago

---

## Quick Start

```bash
# 1. Install dependencies
make setup

# 2. Run full pipeline (loads data + runs all analysis)
make pipeline

# 3. Start the dashboard
make dashboard
```

---

## Dashboard

[View live dashboard](https://pythia-dashboardd.vercel.app) *(replace with deployed dashboard URL)*

---

## How to Run

### Requirements
- Python 3.8+
- `cell-count.csv` in the root directory

### Steps

```bash
# Install dependencies
make setup

# Run pipeline (creates teiko.db and generates all outputs)
make pipeline

# Start dashboard
make dashboard
# Then open http://localhost:8501 in your browser
```

### Running individual scripts

```bash
# Part 1: Initialize database and load data
python load_data.py

# Parts 2-4: Run all analysis
python analyze.py

# Dashboard
streamlit run dashboard.py
```

---

## Database Schema

### Overview

The data is split into three tables to avoid redundancy and support efficient querying.

```
subjects ──────────────────────────────────────────────────────────
  subject    TEXT  PRIMARY KEY   -- unique patient ID
  project    TEXT                -- which clinical trial (prj1, prj2, prj3)
  condition  TEXT                -- melanoma, carcinoma, healthy
  age        INTEGER             -- patient age
  sex        TEXT                -- M or F

samples ───────────────────────────────────────────────────────────
  sample                    TEXT  PRIMARY KEY
  subject                   TEXT  FOREIGN KEY → subjects.subject
  treatment                 TEXT  -- miraclib, phauximab, none
  response                  TEXT  -- yes, no, NULL (healthy patients)
  sample_type               TEXT  -- PBMC or WB
  time_from_treatment_start INTEGER -- 0, 7, or 14 days

cell_counts ───────────────────────────────────────────────────────
  sample      TEXT  PRIMARY KEY  FOREIGN KEY → samples.sample
  b_cell      INTEGER
  cd8_t_cell  INTEGER
  cd4_t_cell  INTEGER
  nk_cell     INTEGER
  monocyte    INTEGER
```

### Why three tables?

**Avoid redundancy.** Each patient has 3 blood draws (time 0, 7, 14). Without normalization, patient-level info like age and sex would repeat 3 times per patient. With 3,500 patients that is 7,000 redundant rows. Separating into a `subjects` table eliminates this.

**Efficient querying.** Analytical queries only join the tables they need. A query about cell counts doesn't need to load subject demographics. This matters at scale.

**Data integrity.** Foreign keys enforce that every sample belongs to a real subject, and every cell count record belongs to a real sample. Orphaned records are impossible.

### How this scales

**Hundreds of projects:** The `project` column in `subjects` handles this cleanly. Adding a separate `projects` table with metadata (start date, indication, sponsor) would be a natural extension without changing the core schema.

**Thousands of samples:** SQLite handles millions of rows without issue. For very large scale (100M+ rows), migrating to PostgreSQL would require no schema changes — just a different connection string.

**Various analytics:**
- Time-series analysis: `time_from_treatment_start` is already indexed by being part of `samples`
- Multi-omics: add new measurement tables (e.g., `protein_counts`, `gene_expression`) that foreign key to `samples`
- New cell populations: add columns to `cell_counts` or create a long-format `measurements` table (sample, population, count) for fully flexible analytics
- Patient cohorts: add a `cohorts` table and a `subject_cohorts` join table

---

## Code Structure

```
├── load_data.py     # Part 1: Creates SQLite database and loads CSV data
├── analyze.py       # Parts 2-4: All analysis functions
├── dashboard.py     # Interactive Streamlit dashboard
├── cell-count.csv   # Input data
├── teiko.db         # Generated SQLite database (created by load_data.py)
├── requirements.txt # Python dependencies
├── Makefile         # Setup, pipeline, and dashboard targets
└── README.md        # This file
```

### load_data.py

Reads `cell-count.csv` and loads it into a normalized 3-table SQLite database. Functions are kept small and single-purpose — `read_csv()`, `create_tables()`, `load_subjects()`, `load_samples()`, `load_cell_counts()`, `verify()`. This makes each step easy to debug and re-run independently.

### analyze.py

Three analysis functions, one per part:

- `compute_cell_frequencies()` — Part 2. Uses `pd.melt()` to transform wide-format cell counts into long format, then calculates percentage per population per sample.
- `compare_responders()` — Part 3. Filters to melanoma + miraclib + PBMC + known response, then runs Mann-Whitney U test for each cell population. Returns filtered data and stats table.
- `analyze_baseline_subset()` — Part 4. SQL-level filtering for baseline samples, then pandas groupby for the four sub-questions.

### dashboard.py

Streamlit app with three tabs — one per analysis part. Each tab shows the data table and relevant visualizations (pie chart, boxplots, bar charts). All data is loaded directly from `teiko.db` at runtime.

### Design philosophy

Code is written to be readable and defensible. Functions are short, named clearly, and commented with the *why* not just the *what*. This matches the requirement to ship readable, maintainable code for high-stakes clinical applications.

---

## Key Findings

**Part 3:** Of the five immune cell populations tested, only **CD4 T cells** showed a statistically significant difference between responders and non-responders (Mann-Whitney U, p=0.0134). Responders had a higher median CD4 T cell frequency (30.22%) compared to non-responders (29.66%), suggesting CD4 T cells may serve as a predictive biomarker for miraclib response in melanoma patients.

**Part 4:** Among 656 baseline melanoma miraclib PBMC samples, the average B cell count for male responders was **10,401.28**.
# teiko-assessment
