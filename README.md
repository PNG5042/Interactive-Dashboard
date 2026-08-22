# Loblaw Bio — Cell Population Analysis

Analysis of immune cell population frequencies from a miraclib clinical trial,
built for Bob Loblaw. Covers database design, sample-level frequency
calculations, a responder vs. non-responder statistical comparison, a
baseline-subset breakdown, and an interactive dashboard tying it together.

## Running it (GitHub Codespaces or local)

```bash
make setup       # installs dependencies from requirements.txt
make pipeline    # builds the database and runs all analysis scripts
make dashboard   # starts the interactive dashboard on port 8501
```

`make pipeline` runs, in order:

1. `load_data.py` — creates `cell_counts.db` and loads `cell-count.csv` into it
2. `initial_analysis.py` — writes `cell_initial_analysis.csv` (Part 2)
3. `statistical_analysis.py` — writes `response_boxplots.png` and `response_stats.csv` (Part 3)
4. `data_subset_analysis.py` — prints the Part 4 breakdown and the bonus B-cell average to stdout

`make dashboard` runs `streamlit run dashboard.py`. In Codespaces, a port
8501 forwarding prompt will pop up automatically — open it to view the
dashboard in your browser. Locally it's at http://localhost:8501.

Dashboard link (hosted): _[add your Streamlit Community Cloud / other hosted URL here]_

## Database schema

Four tables, normalized around the natural entities in the data:

```
projects   (project_id)
subjects   (subject_id, project_id, condition, age, sex, treatment, response)
samples    (sample_id, subject_id, sample_type, time_from_treatment_start)
cell_counts(sample_id, population, count)
```

**Rationale:**

- **One row per subject in `subjects`, not per sample.** Treatment,
  response, condition, age, and sex don't change across a subject's
  samples — storing them once avoids repeating identical values across
  every timepoint and makes updates (e.g. correcting a subject's response
  call) a single-row change instead of a multi-row find-and-replace.
- **`samples` separate from `subjects`.** A subject can contribute multiple
  samples (different timepoints, different sample types), so this is a
  clean one-to-many relationship. `time_from_treatment_start` and
  `sample_type` live here because they vary per sample, not per subject.
- **`cell_counts` in long format** (`sample_id, population, count`) rather
  than wide (one column per population). This is what makes every query in
  this project a simple `GROUP BY population`/`JOIN` instead of hardcoding
  five column names everywhere, and it's the standard shape for this kind
  of data (it's effectively tidy/long-format, which both pandas and SQL
  handle well for aggregation). Adding a sixth population later means
  inserting rows, not migrating the schema.
- **`projects` as its own table** even though it's just an ID today — it's
  a natural place to hang project-level metadata later (PI, institution,
  start date) without touching `subjects`.

**Scaling to hundreds of projects / thousands of samples:**

- The schema itself doesn't need to change — normalization is what makes it
  scale. Wide-format cell counts would require an `ALTER TABLE` every time
  a new population is measured; long-format just gets more rows.
- At real scale, move off a single SQLite file to Postgres (or similar) for
  concurrent writes, and add indexes on `subject_id`, `sample_id`, and
  `population` (SQLite already gets basic ones via the join columns, but a
  bigger dataset would want them explicit and possibly composite).
- For "various types of analytics," this schema already supports most
  slicing (by project, condition, treatment, timepoint, population) via
  joins alone. If specific analytics get expensive at scale (e.g. rolling
  cohort-wide percentiles), a materialized view or a precomputed
  `sample_frequencies` table (basically what `frequency_analysis.py`
  produces) would avoid recomputing percentages on every query.
- Partitioning `cell_counts` by project or by time would help if any single
  query pattern dominates (e.g. always filtering to one project at a time).

## Code structure

```
load_data.py              - Part 1: builds the DB, loads the CSV
initial_analysis.py       - Part 2: per-sample population frequency table
statistical_analysis.py   - Part 3: responder vs non-responder boxplots + stats
data_subset_analysis.py   - Part 4: baseline subset breakdown + bonus B-cell query
dashboard.py               - interactive dashboard (Streamlit)
requirements.txt
Makefile
cell-count.csv             - input data
```

Each analysis script is standalone and re-runnable on its own (`python
initial_analysis.py`, etc.) — they all just read from `cell_counts.db`,
so nothing depends on script execution order except that `load_data.py`
has to run first to create the database. This kept each part of the
assignment isolated and easy to test independently, rather than one long
script doing everything. The dashboard reuses the same DB and re-runs the
same logic live (with an added project filter in Part 2) rather than just
displaying the static outputs, so it stays correct if the underlying data
changes.

**Statistical note (Part 3):** each subject contributes 3 PBMC samples
(different timepoints). Testing at the sample level would pseudoreplicate
— samples from the same subject aren't independent observations — so
`statistical_analysis.py` and the dashboard both average each subject's
samples to one value per population before running the Mann-Whitney U
test, and apply a Benjamini-Hochberg correction across the 5 populations
tested. Result: no population reaches significance after correction
(`cd4_t_cell` is closest, raw p ≈ 0.012, adjusted p ≈ 0.062).
