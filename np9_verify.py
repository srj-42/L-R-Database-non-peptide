#!/usr/bin/env python3
"""Phase NP-9 - verification loop for the chemical half.

Mirrors cur1_verify.py, including its gates, because the chemical half has no
reason to be held to a lower standard than the peptide half - it starts from
weaker sources, so if anything it needs the gates more.

  --emit N     write work/NP_WORK.md: N papers, the pairs each is cited for, and
               the abstract, ordered so one read settles as many pairs as possible
  --load FILE  load adjudicated verdicts, enforcing the gates
  --status     counts by outcome

GATES on a `confirmed` row:
  1. the paper resolves and has text on file
  2. evidence_sentence is VERBATIM in title+abstract
  3. the CHEMICAL is named in that text - by name or by one of its synonyms
  4. the RECEPTOR symbol is named in that text
  5. assay_class is one of the five valid classes

Gate 3 is the one that differs from the peptide loop. A chemical has no gene
symbol, and upstream resources name the same molecule half a dozen ways, so the
match is made against the synonym list carried in chemical_entities - the same
list that did the cross-database joining in np2.

Usage: python3 np9_verify.py --emit 12 | --load work/np_batch1.json | --status
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import unicodedata
import sys

from cur_common import HERE, connect, log, now

WORK = HERE / "work" / "NP_WORK.md"
VALID_OUTCOME = {"confirmed", "pair_absent", "review_only", "abstract_insufficient",
                 "different_entity", "unresolvable", "retracted", "pending"}
VALID_CLASS = {"structural", "direct_binding", "functional_in_vitro",
               "functional_in_vivo", "computational"}


def squash(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "")).strip().lower()


def norm(s: str) -> str:
    """Alphanumerics only - used for the NAME gates, not for the sentence gate.

    Added 2026-10-10. The name gates compared squashed strings, so a paper writing
    'leukotriene C(4)' or 'LTC(4)' never matched the recorded synonyms 'Leukotriene C4'
    and 'LTC4', and 'P2Y14' never matched the symbol 'P2RY14'. The gate then blocked
    confirmations that were correct, which pushes the curator toward editing
    chemical_entities to get past a punctuation mismatch - the wrong fix, because it
    pollutes chemical identity to work around a string comparison. The sentence gate
    still uses squash: a quoted sentence must be verbatim, punctuation included.
    """
    # NFKD folds Unicode subscripts, so a paper's "LTC₄" reads as "LTC4"
    s = unicodedata.normalize("NFKD", s or "")
    return re.sub(r"[^a-z0-9]", "", s.lower())


def chem_names(c, chemical_id: str) -> list[str]:
    r = c.execute("SELECT name, synonyms FROM chemical_entities WHERE chemical_id=?",
                  (chemical_id,)).fetchone()
    if not r:
        return []
    out = [r["name"]]
    try:
        out += [s for s in json.loads(r["synonyms"] or "[]") if s and len(s) > 2]
    except Exception:                                        # noqa: BLE001
        pass
    return out


def emit(n: int) -> int:
    c = connect()
    rows = c.execute("""
        SELECT ne.np_evidence_id, ne.pmid, ne.chem_interaction_id, ne.source_id,
               ch.chemical_id, ch.name chem, ch.hmdb_id, ch.molecular_class,
               e.symbol receptor, ci.receptor_type, ci.n_sources,
               p.title, p.abstract, p.journal, p.year
          FROM np_evidence ne
          JOIN chemical_interactions ci ON ci.chem_interaction_id = ne.chem_interaction_id
          JOIN chemical_entities ch ON ch.chemical_id = ci.chemical_id
          JOIN entities e ON e.entity_id = ci.receptor_entity_id
          JOIN papers p ON p.pmid = ne.pmid
         WHERE ne.verification_outcome = 'pending'
           AND p.abstract IS NOT NULL AND p.abstract != ''
      ORDER BY ne.pmid""").fetchall()
    by_pmid = collections.defaultdict(list)
    for r in rows:
        by_pmid[r["pmid"]].append(r)
    # Most pairs settled per read first.
    order = sorted(by_pmid.items(), key=lambda kv: -len(kv[1]))[:n]

    out = ["# Chemical-half verification work", "",
           f"{len(order)} papers, {sum(len(v) for _, v in order)} pending rows. "
           "Read the abstract; decide each pair.", "",
           "Outcomes: `confirmed` (needs a verbatim sentence + assay_class), `pair_absent` "
           "(tested, does not bind), `abstract_insufficient` (answer is past the abstract), "
           "`different_entity` (wrong molecule, or fails the localisation check), "
           "`review_only`, `unresolvable`.", ""]
    for pmid, rs in order:
        r0 = rs[0]
        out += [f"## PMID {pmid} — {r0['title']}",
                f"*{r0['journal']} {r0['year']}*", "",
                "### Pairs cited to this paper", ""]
        for r in rs:
            syn = chem_names(c, r["chemical_id"])[:6]
            out.append(f"- `{r['np_evidence_id']}` **{r['chem']}** "
                       f"({r['hmdb_id'] or 'no HMDB'}, {r['molecular_class'] or '?'}) "
                       f"→ **{r['receptor']}** ({r['receptor_type'] or 'type?'}) "
                       f"· cited by {r['source_id']} · {r['n_sources']} source(s)")
            out.append(f"    - synonyms: {', '.join(syn)}")
        out += ["", "### Abstract", "", r0["abstract"] or "*none*", "", "---", ""]
    WORK.parent.mkdir(exist_ok=True)
    WORK.write_text("\n".join(out))
    log("np9", f"wrote {WORK} — {len(order)} papers, {sum(len(v) for _, v in order)} rows")
    return 0


def load(path: str) -> int:
    c = connect()
    recs = json.loads(open(path).read())
    counts, gate_log, overwrote = collections.Counter(), [], []
    for r in recs:
        eid = r.get("np_evidence_id")
        row = c.execute("""SELECT ne.np_evidence_id, ne.verification_outcome prior,
                                  ci.chemical_id, e.symbol receptor, p.title, p.abstract,
                                  p.full_text, p.fulltext_status
                             FROM np_evidence ne
                             JOIN chemical_interactions ci
                               ON ci.chem_interaction_id = ne.chem_interaction_id
                             JOIN entities e ON e.entity_id = ci.receptor_entity_id
                             LEFT JOIN papers p ON p.pmid = ne.pmid
                            WHERE ne.np_evidence_id = ?""", (eid,)).fetchone()
        if row is None:
            counts["unknown_id"] += 1
            gate_log.append((eid, "no such np_evidence_id")); continue
        outcome = r.get("verification_outcome")
        if outcome not in VALID_OUTCOME:
            counts["bad_outcome"] += 1
            gate_log.append((eid, f"invalid outcome {outcome!r}")); continue

        if outcome == "confirmed":
            if r.get("evidence_depth") == "full_text":
                # genuine publisher text only (PLAN.md invariants 8-12): same test as cur1 GATE 0
                ft = row["full_text"] or ""
                if row["fulltext_status"] != "pmc_oa" or len(ft) < 20000 or ft.lstrip().startswith("Title:"):
                    counts["gate_not_genuine_full_text"] += 1
                    gate_log.append((eid, "evidence_depth=full_text but no genuine stored full text")); continue
                text = (row["title"] or "") + " " + ft
            else:
                text = (row["title"] or "") + " " + (row["abstract"] or "")
            src = squash(text)
            sent = r.get("evidence_sentence") or ""
            if squash(sent) not in src:
                counts["gate_sentence_not_verbatim"] += 1
                gate_log.append((eid, f"sentence not verbatim: {sent[:60]!r}")); continue
            named = r.get("names_matched") or {}
            nsrc = norm(text)
            cn, rn = norm(named.get("chemical", "")), norm(named.get("receptor", ""))
            valid_chem = {norm(x) for x in chem_names(c, row["chemical_id"])}
            if not cn or cn not in nsrc:
                counts["gate_chemical_not_named"] += 1
                gate_log.append((eid, f"chemical name {named.get('chemical')!r} not in text")); continue
            if cn not in valid_chem:
                counts["gate_chemical_not_a_synonym"] += 1
                gate_log.append((eid, f"{named.get('chemical')!r} is not a recorded name/synonym "
                                      "of this chemical - if it is right, add it in np2")); continue
            if not rn or rn not in nsrc:
                counts["gate_receptor_not_named"] += 1
                gate_log.append((eid, f"receptor name {named.get('receptor')!r} not in text")); continue
            if r.get("assay_class") not in VALID_CLASS:
                counts["gate_bad_assay_class"] += 1
                gate_log.append((eid, f"invalid assay_class {r.get('assay_class')!r}")); continue

        if row["prior"] != "pending":
            overwrote.append((eid, row["prior"], outcome))
        c.execute("""UPDATE np_evidence SET verification_outcome=?, evidence_depth=?,
                        evidence_sentence=?, assay_class=?, effect_direction=?,
                        extraction_method=?, extraction_confidence=?, notes=?, verified_at=?
                      WHERE np_evidence_id=?""",
                  (outcome, r.get("evidence_depth", "abstract_only"),
                   r.get("evidence_sentence"), r.get("assay_class"), r.get("effect_direction"),
                   r.get("extraction_method", "llm_curator_read_abstract"),
                   r.get("extraction_confidence"), r.get("notes"), now(), eid))
        counts[outcome] += 1
    c.commit()
    for k, v in counts.most_common():
        log("np9", f"  {k:32s} {v}")
    for eid, why in gate_log[:20]:
        log("np9", f"  GATE {eid}: {why}")
    if overwrote:
        log("np9", f"  !!! {len(overwrote)} row(s) overwrote a non-pending verdict !!!")
        for eid, a, b in overwrote[:10]:
            log("np9", f"      {eid}: {a} -> {b}")
    return 0


def status() -> int:
    c = connect()
    log("np9", "np_evidence by outcome:")
    for r in c.execute("""SELECT verification_outcome o, count(*) n FROM np_evidence
                          GROUP BY 1 ORDER BY n DESC"""):
        log("np9", f"  {r['o']:24s} {r['n']}")
    n = c.execute("""SELECT count(DISTINCT chem_interaction_id) FROM np_evidence
                      WHERE verification_outcome='confirmed'""").fetchone()[0]
    log("np9", f"chemical pairs with >=1 confirmed paper: {n}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--emit", type=int)
    ap.add_argument("--load")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    sys.exit(emit(a.emit) if a.emit else load(a.load) if a.load else status())
