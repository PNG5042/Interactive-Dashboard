import sqlite3
import pandas as pd

DB_FILE = "cell_counts.db"
OUT_FILE = "cell_initial_analysis.csv"


def compute_frequencies(conn):
    counts = pd.read_sql_query(
        "SELECT sample_id AS sample, population, count FROM cell_counts", conn
    )

    totals = counts.groupby("sample")["count"].sum().rename("total_count")
    counts = counts.join(totals, on="sample")
    counts["percentage"] = counts["count"] / counts["total_count"] * 100

    # keep the requested column order
    counts = counts[["sample", "total_count", "population", "count", "percentage"]]
    counts = counts.sort_values(["sample", "population"]).reset_index(drop=True)
    return counts


if __name__ == "__main__":
    conn = sqlite3.connect(DB_FILE)
    freq_table = compute_frequencies(conn)
    conn.close()

    freq_table.to_csv(OUT_FILE, index=False)

    print(f"wrote {len(freq_table)} rows to {OUT_FILE}")
    print()
    print(freq_table.head(10).to_string(index=False))