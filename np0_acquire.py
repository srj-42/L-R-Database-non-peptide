#!/usr/bin/env python3
"""Phase NP-0 - acquire and PIN the non-peptide sources.

Same contract as cur0a_acquire.py: download if absent, then record version, URL,
SHA-256, row count, download date and licence in `sources`. Adds one thing the
peptide phase did not need - a `shippable` flag, because MEBOCOST's licence
forbids redistribution and that has to be enforceable at export time.

Sources already pinned by cur0a (GtoPdb, CellChatDB, CellPhoneDB) are not
re-downloaded here; np0 only asserts they are present and sets their flag.

See PHASE_NONPEPTIDE.md for the decisions this implements.
"""
from __future__ import annotations
import json, sqlite3, sys
from pathlib import Path

from cur_common import (HERE, RAW, connect, init_schema, log, now, read_csv,
                        sha256, fetch)

SRC_NP = json.loads((HERE / "sources_np.json").read_text())


def migrate(c: sqlite3.Connection) -> None:
    """Add sources.shippable if it is not there yet. Idempotent."""
    cols = {r["name"] for r in c.execute("PRAGMA table_info(sources)")}
    if "shippable" not in cols:
        c.execute("ALTER TABLE sources ADD COLUMN shippable INTEGER DEFAULT 1")
        c.commit()
        log("np0", "migrated: sources.shippable added")


def row_count(p: Path, delim: str | None) -> int | None:
    """Row count for the formats this phase actually uses."""
    if p.suffix == ".db":
        try:
            with sqlite3.connect(p) as s:
                return s.execute("SELECT count(*) FROM edges").fetchone()[0]
        except Exception:                                  # noqa: BLE001
            return None
    if p.suffix in (".csv", ".tsv", ".txt") and delim:
        try:
            return len(read_csv(p, delim))
        except Exception:                                  # noqa: BLE001
            return None
    return None


def pin(c: sqlite3.Connection, sid: str, meta: dict, dest: Path) -> None:
    digest = sha256(dest)
    expected = meta.get("sha256_expected")
    if expected and digest != expected:
        raise RuntimeError(
            f"{sid}: sha256 mismatch\n  expected {expected}\n  got      {digest}\n"
            "The upstream release asset changed. Do not proceed silently - the "
            "counts in PHASE_NONPEPTIDE.md were measured against the pinned file.")
    n = row_count(dest, meta.get("delim"))
    c.execute(
        "INSERT OR REPLACE INTO sources"
        "(source_id,name,version,url,local_file,sha256,n_rows,downloaded_at,licence,role,shippable)"
        " VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        (sid, meta["name"], None, meta["url"], str(dest.relative_to(HERE.parent)),
         digest, n, now(), meta["licence"], meta["role"], meta["shippable"]))
    flag = "" if meta["shippable"] else "  [NOT SHIPPABLE]"
    log("np0", f"{sid:30s} rows={n if n is not None else '-':>7} sha={digest[:12]}{flag}")


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    c = connect()
    init_schema(c)
    migrate(c)

    # Sources the peptide phase already pinned. Assert, do not re-fetch.
    missing = []
    for sid, meta in SRC_NP["reuse"].items():
        row = c.execute("SELECT source_id FROM sources WHERE source_id=?", (sid,)).fetchone()
        if row is None:
            missing.append(sid)
            continue
        c.execute("UPDATE sources SET shippable=? WHERE source_id=?", (meta["shippable"], sid))
        log("np0", f"{sid:30s} reused from sources.json")
    if missing:
        log("np0", f"MISSING from `sources`: {', '.join(missing)}")
        log("np0", "run cur0a_acquire.py first - np0 does not re-pin peptide-phase sources")
        return 1

    for sid, meta in SRC_NP["sources"].items():
        dest = (HERE / meta["local_file"]).resolve()
        if not dest.exists():
            log("np0", f"downloading {sid} ...")
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(fetch(meta["url"]))
        pin(c, sid, meta, dest)

    c.commit()

    n_ship = c.execute("SELECT count(*) FROM sources WHERE shippable=0").fetchone()[0]
    log("np0", f"pinned {len(SRC_NP['sources'])} new + {len(SRC_NP['reuse'])} reused sources")
    log("np0", f"{n_ship} source(s) marked NOT shippable - np6_export.py must honour this")
    return 0


if __name__ == "__main__":
    sys.exit(main())
