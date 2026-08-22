import sqlite3
import pandas as pd

DB_FILE = "cell_counts.db"


def get_baseline_subset(conn):
    query = """
        SELECT
            sm.sample_id,
            sm.subject_id,
            s.project_id,
            s.response,
            s.sex
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        WHERE s.condition = 'melanoma'
          AND s.treatment = 'miraclib'
          AND sm.sample_type = 'PBMC'
          AND sm.time_from_treatment_start = 0
    """
    return pd.read_sql_query(query, conn)


if __name__ == "__main__":
    conn = sqlite3.connect(DB_FILE)
    subset = get_baseline_subset(conn)
    conn.close()

    print(f"baseline melanoma/miraclib/PBMC samples: {len(subset)}")
    print(f"unique subjects: {subset['subject_id'].nunique()}")
    print()

    # one row per sample here, but since this is baseline (t=0) each
    # subject only shows up once, so sample counts == subject counts
    print("samples per project:")
    print(subset["project_id"].value_counts().to_string())
    print()

    subjects = subset.drop_duplicates("subject_id")

    print("responders vs non-responders:")
    print(subjects["response"].value_counts().to_string())
    print()

    print("male vs female:")
    print(subjects["sex"].value_counts().to_string())
    print()

    print("-" * 60)
    avg_query = """
        SELECT AVG(cc.count) AS avg_b_cell
        FROM samples sm
        JOIN subjects s ON s.subject_id = sm.subject_id
        JOIN cell_counts cc ON cc.sample_id = sm.sample_id
        WHERE s.condition = 'melanoma'
          AND s.sex = 'M'
          AND s.response = 'yes'
          AND sm.time_from_treatment_start = 0
          AND cc.population = 'b_cell'
    """
    conn = sqlite3.connect(DB_FILE)
    avg_b_cell = pd.read_sql_query(avg_query, conn).iloc[0]["avg_b_cell"]
    conn.close()
    print(f"avg B cell count, melanoma male responders, all sample/treatment types, t=0: {avg_b_cell:.2f}")