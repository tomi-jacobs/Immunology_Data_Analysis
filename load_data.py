"""
load_data.py
Tomi Jacobs
May 2026

Purpose:
    Read cell-count.csv and load it into a SQLite database
    with 3 tables: subjects, samples, cell_counts
"""

import sqlite3
import pandas as pd
import os


# -- File names --
CSV_FILE = "cell-count.csv"
DB_FILE = "teiko.db"


# Step 1: Read the CSV file into a DataFrame
def read_csv():
    """Read the CSV file"""
    data = pd.read_csv(CSV_FILE)
    print(f"Read {len(data)} rows from {CSV_FILE}")
    return data


# Step 2: Create the 3 tables in SQLite
def create_tables(conn):
    """
    Create 3 tables: subjects, samples, cell_counts

    Why 3 tables and not 1?
    Each patient has 3 blood draws (time 0, 7, 14).
    If we use 1 table, patient info like age and sex
    repeats 3 times. Splitting into 3 tables avoids
    this repetition. This is called normalization.

    subjects    --> one row per patient
    samples     --> one row per blood draw
    cell_counts --> one row per blood draw with cell measurements
    """

    cursor = conn.cursor()

    # Drop old tables so the script can be re-run safely
    cursor.execute("DROP TABLE IF EXISTS cell_counts")
    cursor.execute("DROP TABLE IF EXISTS samples")
    cursor.execute("DROP TABLE IF EXISTS subjects")

    # Table 1: subjects
    # Stores patient-level info that does not change over time
    cursor.execute("""
        CREATE TABLE subjects (
            subject   TEXT PRIMARY KEY,
            project   TEXT,
            condition TEXT,
            age       INTEGER,
            sex       TEXT
        )
    """)

    # Table 2: samples
    # Stores info specific to each blood draw
    # subject links back to the subjects table (foreign key)
    cursor.execute("""
        CREATE TABLE samples (
            sample                    TEXT PRIMARY KEY,
            subject                   TEXT,
            treatment                 TEXT,
            response                  TEXT,
            sample_type               TEXT,
            time_from_treatment_start INTEGER,
            FOREIGN KEY (subject) REFERENCES subjects(subject)
        )
    """)

    # Table 3: cell_counts
    # Stores the 5 immune cell measurements per blood draw
    # sample links back to the samples table (foreign key)
    cursor.execute("""
        CREATE TABLE cell_counts (
            sample     TEXT PRIMARY KEY,
            b_cell     INTEGER,
            cd8_t_cell INTEGER,
            cd4_t_cell INTEGER,
            nk_cell    INTEGER,
            monocyte   INTEGER,
            FOREIGN KEY (sample) REFERENCES samples(sample)
        )
    """)

    conn.commit()
    print("Tables created successfully")


# Step 3: Load subjects into the subjects table
def load_subjects(conn, data):
    """
    Insert one row per unique patient.
    Each patient appears 3 times in the CSV (one per time point).
    drop_duplicates() keeps only the first occurrence of each patient.
    """

    # Select only the columns needed for this table
    subjects = data[["subject", "project", "condition", "age", "sex"]]

    # Remove duplicate patients
    subjects = subjects.drop_duplicates(subset="subject")

    # Insert into database
    subjects.to_sql("subjects", conn, if_exists="append", index=False)

    print(f"Loaded {len(subjects)} subjects")


# Step 4: Load samples into the samples table
def load_samples(conn, data):
    """
    Insert one row per blood draw.
    Healthy patients have no response value - stored as NULL in SQLite.
    """

    # Select only the columns needed for this table
    samples = data[["sample", "subject", "treatment", "response",
                    "sample_type", "time_from_treatment_start"]].copy()

    # NaN in pandas becomes NULL in SQLite automatically
    samples.to_sql("samples", conn, if_exists="append", index=False)

    print(f"Loaded {len(samples)} samples")


# Step 5: Load cell counts into the cell_counts table
def load_cell_counts(conn, data):
    """
    Insert one row per blood draw with all 5 cell measurements.
    """

    # Select only the columns needed for this table
    cell_counts = data[["sample", "b_cell", "cd8_t_cell",
                         "cd4_t_cell", "nk_cell", "monocyte"]]

    cell_counts.to_sql("cell_counts", conn, if_exists="append", index=False)

    print(f"Loaded {len(cell_counts)} cell count records")


# Step 6: Verify everything loaded correctly
def verify(conn):
    """Print row counts for each table"""

    cursor = conn.cursor()

    print("\nVerification:")
    for table in ["subjects", "samples", "cell_counts"]:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        count = cursor.fetchone()[0]
        print(f"  {table}: {count} rows")


# -- Main: run all steps in order --
def main():

    # Remove old database if it exists
    if os.path.exists(DB_FILE):
        os.remove(DB_FILE)
        print(f"Removed old {DB_FILE}")

    # Read the CSV file
    data = read_csv()

    # Connect to SQLite database
    conn = sqlite3.connect(DB_FILE)

    # Create the 3 tables
    create_tables(conn)

    # Load data into each table
    load_subjects(conn, data)
    load_samples(conn, data)
    load_cell_counts(conn, data)

    # Verify row counts
    verify(conn)

    # Close connection
    conn.close()

    print(f"\nDone! Database saved as {DB_FILE}")


# Run main when script is executed directly
if __name__ == "__main__":
    main()
