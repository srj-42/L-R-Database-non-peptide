#!/usr/bin/env python3
"""Merge name-keyed duplicate chemical_entities into their identified twins.

The gtopdb_interactions ingest created chemical rows keyed by name only - no
hmdb_id, pubchem_cid or chebi_id - for compounds already present with an
identifier from the other resources. 'ATP' (CHEM00630, no ids, gtopdb only)
duplicates 'Adenosine triphosphate' (CHEM00149, HMDB0000538), and so on. No
identifier-based check can see this, because the duplicates have no identifiers;
PLAN.md's note that "the remaining 99 are name-keyed only and are exactly the
ones most likely to be confused with a near neighbour" is this failure mode.

Effect of leaving them: 83 human chemical_interactions rows are duplicates of
another pair, so every downstream pass that treats pairs as independent -
scoring, export, this phase's literature search - does the work twice and
double-counts the result.

MERGES is hand-reviewed, not generated. Synonym matching proposed the
candidates, but it CANNOT be trusted to apply them: CHEM00118
(UDP-glucose, HMDB0000286) carries the bare synonym "UDP", which wrongly
matches CHEM00687 (UDP). That one is mapped to CHEM00122 (Uridine
5'-diphosphate, HMDB0000295) below, and ambiguous candidates are refused.

Three merges fold a stereo-specific name into an unspecified parent whose HMDB
entry IS that stereoisomer (12S-HETE/12-HETE HMDB0006111, 8S-HETE/8-HETE
HMDB0004679, 24(S)-hydroxycholesterol/24-Hydroxycholesterol HMDB0001419). They
are listed separately as STEREO_REVIEWED so a curator can see what was assumed.

What a merge does, per duplicate:
  * chemical_interactions: a dup pair whose (keep, receptor, species) row already
    exists is folded into it - sources unioned, n_sources recomputed,
    has_citation OR-ed, n_pmids maxed - and the dup row deleted; otherwise the
    row is simply repointed to the kept chemical_id, preserving its
    chem_interaction_id and therefore its np_evidence.
  * np_evidence: repointed to the surviving chem_interaction_id. A row that
    would collide on UNIQUE(chem_interaction_id, pmid) is dropped, because the
    surviving row already cites that paper for that pair.
  * chem_source_records, metabolite_enzymes: repointed; metabolite_enzymes
    collisions on UNIQUE(chemical_id, gene_symbol, role, source_id) are dropped.
  * chemical_entities: the duplicate row is deleted last.

Nothing here touches papers.* or any evidence_sentence (invariant 8).

Usage: python3 np_dedup_chemicals.py [--commit]
"""
from __future__ import annotations

import json
import sys

from cur_common import connect, log

# duplicate (name-keyed, no identifier) -> keep (carries the identifier)
MERGES = {
    "CHEM00600": "CHEM00278",  # 20-hydroxy-LTB4            -> 20-Hydroxy-leukotriene B4
    "CHEM00608": "CHEM00268",  # 5a-pregnan-3a-ol-20-one    -> Allopregnanolone
    "CHEM00609": "CHEM00108",  # 5-hydroxytryptamine        -> Serotonin
    "CHEM00613": "CHEM00412",  # 5-oxo-ETE                  -> 5-KETE
    "CHEM00617": "CHEM00276",  # 7a-hydroxycholesterol      -> 7alpha-Hydroxycholesterol
    "CHEM00620": "CHEM00389",  # 7b-hydroxycholesterol      -> 7beta-Hydroxycholesterol
    "CHEM00626": "CHEM00254",  # a-linolenic acid           -> alpha-Linolenic acid
    "CHEM00627": "CHEM00314",  # alitretinoin               -> 9-cis-Retinoic acid
    "CHEM00630": "CHEM00149",  # ATP                        -> Adenosine triphosphate
    "CHEM00631": "CHEM00024",  # b-alanine                  -> beta-Alanine
    "CHEM00633": "CHEM00422",  # b-phenylethylamine         -> Phenylethylamine
    "CHEM00639": "CHEM00112",  # dihydrosphingosine         -> Sphinganine
    "CHEM00644": "CHEM00212",  # farnesyl diphosphate       -> Farnesyl pyrophosphate
    "CHEM00649": "CHEM00333",  # H2O2                       -> Hydrogen peroxide
    "CHEM00653": "CHEM00206",  # L-citrulline               -> Citrulline
    "CHEM00654": "CHEM00077",  # levodopa                   -> L-Dopa
    "CHEM00655": "CHEM00087",  # L-ornithine                -> Ornithine
    "CHEM00661": "CHEM00225",  # methylglyoxal              -> Pyruvaldehyde
    "CHEM00676": "CHEM00187",  # pregnenolone sulphate      -> Pregnenolone sulfate
    "CHEM00680": "CHEM00127",  # retinol                    -> Vitamin A (HMDB0000305 = retinol)
    "CHEM00683": "CHEM00102",  # T4                         -> Thyroxine
    "CHEM00685": "CHEM00212",  # trans,trans-farnesyl diphosphate -> Farnesyl pyrophosphate
    "CHEM00686": "CHEM00284",  # tretinoin                  -> all-trans-Retinoic acid
    "CHEM00687": "CHEM00122",  # UDP -> Uridine 5'-diphosphate   NOT CHEM00118 (see docstring)
    "CHEM00688": "CHEM00125",  # UDP-galactose              -> Uridine diphosphategalactose
    "CHEM00689": "CHEM00118",  # UDP-glucose                -> Uridine diphosphate glucose
    "CHEM00690": "CHEM00117",  # UTP                        -> Uridine triphosphate
}

# Round 2, 2026-10-10. Found while auditing upstream citations: np9's chemical gate
# refused correct confirmations because the paper's spelling was not a recorded name, and
# the collision guard in np2c_add_abbreviations.py then showed WHY - the spelling was
# already a name, of a second entity holding the same molecule. Detection by name alone
# missed these, because 'Geranyl-PP' and 'geranyl diphosphate' share no normalised form;
# they surfaced only when a synonym was proposed and found to be taken.
MERGES_R2 = {
    "CHEM00641": "CHEM00222",  # dimethylallyl diphosphate  -> Dimethylallylpyrophosphate
    "CHEM00646": "CHEM00240",  # geranyl diphosphate        -> Geranyl-PP
    "CHEM00652": "CHEM00248",  # isopentenyl diphosphate    -> Isopentenyl pyrophosphate
    "CHEM00647": "CHEM00362",  # geranylgeranyl diphosphate -> Geranylgeranyl-PP
    "CHEM00677": "CHEM00096",  # propanoic acid             -> Propionic acid
    "CHEM00658": "CHEM00220",  # LTB4                       -> Leukotriene B4
    "CHEM00659": "CHEM00230",  # LTC4                       -> Leukotriene C4
    "CHEM00672": "CHEM00223",  # PGF2a                      -> Prostaglandin F2a
    "CHEM00673": "CHEM00245",  # PGI2                       -> Prostaglandin I2
}

# Round 3, 2026-10-10. Same finding as round 2 - the audit kept stalling on spellings
# ('PGJ(2)', 'Ca(2+)', 'IP3') that belonged to a second, name-keyed entity. Only pairs
# where both records are certainly one molecule are listed. Deliberately NOT merged:
# 15S-HETE / 5S-HETE (the identified 15-HETE / 5-HETE do not state a stereocentre),
# Cu2+ (oxidation state not fixed by 'Copper'), the aspirin-triggered lipoxin A4 and
# resolvin D1 (17R / 15-epi stereoisomers of the identified compounds), 7b,25-diOHC.
MERGES_R3 = {
    "CHEM00598": "CHEM00060",  # 17b-estradiol               -> Estradiol
    "CHEM00591": "CHEM00288",  # 1,25-dihydroxyvitamin D3    -> Calcitriol
    "CHEM00599": "CHEM00386",  # 20-HETE                     -> 20-Hydroxyeicosatetraenoic acid
    "CHEM00624": "CHEM00226",  # ADP ribose                  -> Adenosine diphosphate ribose
    "CHEM00623": "CHEM00019",  # adenosine 5'-monophosphate  -> Adenosine monophosphate
    "CHEM00634": "CHEM00138",  # Ca2+                        -> Calcium
    "CHEM00662": "CHEM00150",  # Mg2+                        -> Magnesium
    "CHEM00663": "CHEM00157",  # Na+                         -> Sodium
    "CHEM00691": "CHEM00523",  # Zn2+                        -> Zinc
    "CHEM00651": "CHEM00277",  # IP3                         -> myo-Inositol 1,4,5-trisphosphate
    "CHEM00665": "CHEM00219",  # NAAG                        -> N-Acetylaspartylglutamic acid
    "CHEM00669": "CHEM00583",  # PAF                         -> Platelet-activating factor
    "CHEM00670": "CHEM00266",  # PGE1                        -> Prostaglandin E1
    "CHEM00671": "CHEM00317",  # PGF1a                       -> Prostaglandin F1a
    "CHEM00674": "CHEM00319",  # PGJ2                        -> Prostaglandin J2
    "CHEM00667": "CHEM00303",  # N-oleoylethanolamide        -> Oleoylethanolamide
    "CHEM00668": "CHEM00304",  # N-palmitoylethanolamine     -> Palmitoylethanolamide
    "CHEM00640": "CHEM00253",  # dihydrosphingosine 1-phosphate -> Sphinganine 1-phosphate
    "CHEM00648": "CHEM00580",  # H+                          -> Hydrogen Ion
    "CHEM00625": "CHEM00086",  # a-ketoglutaric acid         -> Oxoglutaric acid
    "CHEM00632": "CHEM00001",  # b-D-hydroxybutyric acid -> 3-Hydroxybutyric acid ((3R) = D)
    "CHEM00657": "CHEM00114",  # L-serine-O-phosphate        -> Phosphoserine (O-phospho-L-serine)
    "CHEM00596": "CHEM00378",  # 15-deoxy-d-12,14-prostaglandin J2 -> 15-Deoxy-d-12,14-PGJ2
    "CHEM00678": "CHEM00221",  # (R)-5-diphosphomevalonate  -> (R)-5-Diphosphomevalonic acid
    "CHEM00593": "CHEM00423",  # 12-hydroxyheptadecatrienoic acid -> 12S-HHT (12-HHT is the natural 12(S) isomer)
}

# Same merge, but the kept row's name does not state the stereochemistry the
# duplicate's name does. Its HMDB entry is that stereoisomer, so these are
# folded in too - flagged separately for a curator read.
STEREO_REVIEWED = {
    "CHEM00595": "CHEM00388",  # 12S-HETE                  -> 12-HETE (HMDB0006111)
    "CHEM00621": "CHEM00367",  # 8S-HETE                   -> 8-HETE  (HMDB0004679)
    "CHEM00604": "CHEM00261",  # 24(S)-hydroxycholesterol  -> 24-Hydroxycholesterol (HMDB0001419)
}


def jload(s):
    try:
        v = json.loads(s or "[]")
        return v if isinstance(v, list) else [v]
    except Exception:
        return [s] if s else []


def main() -> int:
    commit = "--commit" in sys.argv
    c = connect()
    plan = dict(MERGES) | dict(MERGES_R2) | dict(MERGES_R3) | dict(STEREO_REVIEWED)

    # Idempotent: a duplicate that no longer exists was merged by an earlier run, so it
    # is dropped from the plan rather than aborting it. A missing TARGET is still fatal,
    # because that would mean merging into something that is not there.
    for dup in [d for d in plan if not c.execute(
            "SELECT 1 FROM chemical_entities WHERE chemical_id=?", (d,)).fetchone()]:
        log("dedup", f"  {dup} already merged by an earlier run, skipping")
        plan.pop(dup)
    for dup, keep in plan.items():
        if not c.execute("SELECT 1 FROM chemical_entities WHERE chemical_id=?", (keep,)).fetchone():
            log("dedup", f"ABORT: merge target {keep} does not exist")
            return 2
        if keep in plan:
            log("dedup", f"ABORT: target {keep} is itself scheduled for merge")
            return 2
        row = c.execute("SELECT hmdb_id, pubchem_cid, chebi_id FROM chemical_entities "
                        "WHERE chemical_id=?", (keep,)).fetchone()
        if not (row["hmdb_id"] or row["pubchem_cid"] or row["chebi_id"]):
            log("dedup", f"ABORT: target {keep} has no identifier - not a valid merge target")
            return 2

    st = dict(folded=0, repointed=0, ev_moved=0, ev_dropped=0, ev_source_kept=0, src_moved=0,
              enz_moved=0, enz_dropped=0, chem_deleted=0)

    for dup, keep in plan.items():
        for ci in c.execute("SELECT * FROM chemical_interactions WHERE chemical_id=?", (dup,)).fetchall():
            tgt = c.execute("""SELECT * FROM chemical_interactions
                                WHERE chemical_id=? AND receptor_entity_id=? AND species=?""",
                            (keep, ci["receptor_entity_id"], ci["species"])).fetchone()
            if tgt:
                srcs = sorted(set(jload(tgt["sources"])) | set(jload(ci["sources"])))
                c.execute("""UPDATE chemical_interactions
                                SET sources=?, n_sources=?, has_citation=?, n_pmids=?
                              WHERE chem_interaction_id=?""",
                          (json.dumps(srcs), len(srcs),
                           max(tgt["has_citation"], ci["has_citation"]),
                           max(tgt["n_pmids"], ci["n_pmids"]), tgt["chem_interaction_id"]))
                for ev in c.execute("SELECT np_evidence_id, pmid, source_id FROM np_evidence "
                                    "WHERE chem_interaction_id=?", (ci["chem_interaction_id"],)).fetchall():
                    clash = c.execute("SELECT np_evidence_id, source_id, notes FROM np_evidence "
                                      "WHERE chem_interaction_id=? AND pmid=?",
                                      (tgt["chem_interaction_id"], ev["pmid"])).fetchone()
                    if clash:
                        # Same pair, same paper: the survivor already covers it. But 144 of
                        # these carry a DIFFERENT source_id (gtopdb_interactions vs the
                        # resource that produced the survivor), and that corroboration is
                        # what n_sources is built from. Record it on the survivor as a note
                        # PREFIX before deleting - a prefix keeps any trailing
                        # "Prompt answer: {json}" parseable.
                        if (ev["source_id"] or "") != (clash["source_id"] or ""):
                            mark = (f"[dedup {dup}->{keep}] this pair/paper was also cited by "
                                    f"{ev['source_id']} (duplicate evidence row "
                                    f"{ev['np_evidence_id']} removed). ")
                            if mark not in (clash["notes"] or ""):
                                c.execute("UPDATE np_evidence SET notes=? WHERE np_evidence_id=?",
                                          (mark + (clash["notes"] or ""), clash["np_evidence_id"]))
                                st["ev_source_kept"] += 1
                        c.execute("DELETE FROM np_evidence WHERE np_evidence_id=?", (ev["np_evidence_id"],))
                        st["ev_dropped"] += 1
                    else:
                        c.execute("UPDATE np_evidence SET chem_interaction_id=? WHERE np_evidence_id=?",
                                  (tgt["chem_interaction_id"], ev["np_evidence_id"]))
                        st["ev_moved"] += 1
                c.execute("UPDATE chem_source_records SET chem_interaction_id=?, chemical_id=? "
                          "WHERE chem_interaction_id=?",
                          (tgt["chem_interaction_id"], keep, ci["chem_interaction_id"]))
                st["src_moved"] += c.execute("SELECT changes()").fetchone()[0]
                c.execute("DELETE FROM chemical_interactions WHERE chem_interaction_id=?",
                          (ci["chem_interaction_id"],))
                st["folded"] += 1
            else:
                c.execute("UPDATE chemical_interactions SET chemical_id=? WHERE chem_interaction_id=?",
                          (keep, ci["chem_interaction_id"]))
                c.execute("UPDATE chem_source_records SET chemical_id=? WHERE chem_interaction_id=?",
                          (keep, ci["chem_interaction_id"]))
                st["repointed"] += 1

        # any chem_source_records still pointing at the dup chemical
        c.execute("UPDATE chem_source_records SET chemical_id=? WHERE chemical_id=?", (keep, dup))

        for me in c.execute("SELECT id, gene_symbol, role, source_id FROM metabolite_enzymes "
                            "WHERE chemical_id=?", (dup,)).fetchall():
            clash = c.execute("""SELECT 1 FROM metabolite_enzymes
                                  WHERE chemical_id=? AND gene_symbol=? AND role=? AND source_id=?""",
                              (keep, me["gene_symbol"], me["role"], me["source_id"])).fetchone()
            if clash:
                c.execute("DELETE FROM metabolite_enzymes WHERE id=?", (me["id"],))
                st["enz_dropped"] += 1
            else:
                c.execute("UPDATE metabolite_enzymes SET chemical_id=? WHERE id=?", (keep, me["id"]))
                st["enz_moved"] += 1

        c.execute("DELETE FROM chemical_entities WHERE chemical_id=?", (dup,))
        st["chem_deleted"] += 1

    for k, v in st.items():
        log("dedup", f"  {k:14s} {v}")
    left = c.execute("SELECT COUNT(*) FROM chemical_interactions WHERE species='human'").fetchone()[0]
    log("dedup", f"human chemical_interactions now {left}")
    (c.commit if commit else c.rollback)()
    log("dedup", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
