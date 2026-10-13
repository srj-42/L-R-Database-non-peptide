"""Shared guards for the Phase NP build. Standard library only.

Written after a silent-failure incident on 2026-09-24: `SPECIES = "human"` was
compared against GtoPdb's `Target Species` value `"Human"`, so every GtoPdb row
was skipped. It raised no error, wrote no dropped-row entry and printed no zero -
the source simply did not appear in the output, and the build looked healthy.
GtoPdb is the only source in the union with expert curation and quantitative
affinities, so shipping without it would have been a serious quiet loss.

Two defences live here.

1. CASE-INSENSITIVE COMPARISON of the enum-like strings sources use. Their
   capitalisation is not a stable interface: 'Human'/'human', 'true'/'TRUE',
   'FALSE'/'False' all appear across the seven files, and a release can change
   one without notice.

2. CONTRIBUTION GUARDS. Every phase records how many rows each source
   contributed, and compares against the previous successful run. A source that
   contributes nothing, or collapses relative to last time, stops the build.
   Comparing against history rather than a hardcoded floor means the guard keeps
   working as the sources are updated.
"""
from __future__ import annotations

import sqlite3

from cur_common import log, now

# Fraction of the previous run below which a source is considered to have
# collapsed. 0.5 is deliberately loose - it is there to catch a filter that
# broke, not to police normal upstream churn.
COLLAPSE_RATIO = 0.5

CONTRIB_DDL = """
CREATE TABLE IF NOT EXISTS np_run_contributions (
    phase     TEXT NOT NULL,
    source_id TEXT NOT NULL,
    n_rows    INTEGER NOT NULL,
    run_at    TEXT NOT NULL,
    PRIMARY KEY (phase, source_id)
);
"""


def eqi(a, b) -> bool:
    """Case-insensitive equality for source-provided enum strings."""
    return str(a or "").strip().lower() == str(b or "").strip().lower()


def in_ci(value, options) -> bool:
    """Case-insensitive membership for source-provided enum strings."""
    return str(value or "").strip().lower() in {str(o).strip().lower() for o in options}


def guard_contributions(c: sqlite3.Connection, phase: str,
                        counts: dict[str, int], expected: set[str],
                        commit: bool) -> list[str]:
    """Compare this run's per-source contribution against the last one.

    Returns a list of problem descriptions; empty means the build may proceed.
    Records the new counts only when `commit` is set and there are no problems,
    so a broken run never becomes the baseline a later run is judged against.
    """
    c.executescript(CONTRIB_DDL)
    previous = {r[0]: r[1] for r in c.execute(
        "SELECT source_id, n_rows FROM np_run_contributions WHERE phase=?", (phase,))}

    problems: list[str] = []
    for sid in sorted(expected):
        n = counts.get(sid, 0)
        prev = previous.get(sid)
        if n == 0:
            problems.append(
                f"{sid}: contributed 0 rows"
                + (f" (was {prev} last run)" if prev else " (no previous run to compare)")
                + " - a source contributing nothing is a bug, not an empty source")
        elif prev and n < prev * COLLAPSE_RATIO:
            problems.append(
                f"{sid}: contributed {n} rows, down from {prev} "
                f"({n / prev:.0%} of the previous run) - check the filters before proceeding")

    unexpected = sorted(set(counts) - expected)
    for sid in unexpected:
        log(phase, f"NOTE: {sid} contributed {counts[sid]} rows but is not in the "
                   "expected set - add it there if that is intended")

    if problems:
        for p in problems:
            log(phase, f"ABORT: {p}")
        return problems

    for sid in sorted(expected):
        n = counts.get(sid, 0)
        prev = previous.get(sid)
        delta = "" if prev is None else (f"  (was {prev})" if prev != n else "  (unchanged)")
        log(phase, f"  contributed {sid:28s} {n:6d}{delta}")
    if commit:
        c.executemany(
            "INSERT OR REPLACE INTO np_run_contributions (phase,source_id,n_rows,run_at)"
            " VALUES (?,?,?,?)",
            [(phase, sid, counts.get(sid, 0), now()) for sid in expected])
    return []


# Every source Phase NP reads. A phase passes the subset it is meant to touch.
NP_SOURCES = {
    "mrclinkdb_human_mr",
    "metaligand_pairs_human",
    "metalinksdb",
    "mebocost_human",
    "cellchatdb_human",
    "cellphonedb_interactions",
    "gtopdb_interactions",
}
