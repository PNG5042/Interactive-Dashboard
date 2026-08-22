import csv
import os
import sqlite3

CSV_FILE = "cell-count.csv"
DB_FILE = "cell_counts.db"

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
CREATE TABLE projects (
    project_id TEXT PRIMARY KEY
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT REFERENCES projects(project_id),
    condition TEXT,
    age INTEGER,
    sex TEXT,
    treatment TEXT,
    response TEXT
);

CREATE TABLE samples (
    sample_id TEXT PRIMARY KEY,
    subject_id TEXT REFERENCES subjects(subject_id),
    sample_type TEXT,
    time_from_treatment_start REAL
);

CREATE TABLE cell_counts (
    sample_id TEXT REFERENCES samples(sample_id),
    population TEXT,
    count INTEGER,
    PRIMARY KEY (sample_id, population)
);
"""


def build_db():
    if os.path.exists(DB_FILE):
        os.remove(DB_FILE)

    conn = sqlite3.connect(DB_FILE)
    cur = conn.cursor()
    cur.executescript(SCHEMA)

    seen_projects = set()
    seen_subjects = set()

    with open(CSV_FILE, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_id = row["project"]
            subject_id = row["subject"]
            sample_id = row["sample"]

            if project_id not in seen_projects:
                cur.execute("INSERT INTO projects VALUES (?)", (project_id,))
                seen_projects.add(project_id)

            if subject_id not in seen_subjects:
                cur.execute(
                    "INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        subject_id,
                        project_id,
                        row["condition"],
                        int(row["age"]),
                        row["sex"],
                        row["treatment"],
                        row["response"] if row["response"] else None,
                    ),
                )
                seen_subjects.add(subject_id)

            cur.execute(
                "INSERT INTO samples VALUES (?, ?, ?, ?)",
                (sample_id, subject_id, row["sample_type"], float(row["time_from_treatment_start"])),
            )

            for pop in POPULATIONS:
                cur.execute(
                    "INSERT INTO cell_counts VALUES (?, ?, ?)",
                    (sample_id, pop, int(row[pop])),
                )

    conn.commit()

    # quick sanity check on the way out
    n_samples = cur.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    n_subjects = cur.execute("SELECT COUNT(*) FROM subjects").fetchone()[0]
    n_counts = cur.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0]
    conn.close()

    print(f"done -> {DB_FILE}")
    print(f"  {len(seen_projects)} projects, {n_subjects} subjects, {n_samples} samples, {n_counts} cell count rows")


if __name__ == "__main__":
    build_db()