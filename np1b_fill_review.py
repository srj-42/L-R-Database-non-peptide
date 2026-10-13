#!/usr/bin/env python3
"""Phase NP-1b - first pass at the np1 review queue.

Fills `manual_chem_name`, `manual_hmdb_id` and `manual_enzyme_genes` in
np_composite_map.tsv for the 37 rows np1 could not resolve.

Two rules this script follows, and they matter:

  1. The ABBREVIATION EXPANSION is a human judgement, written out below with a
     rationale per entry. It is reviewable text, not a lookup.
  2. The HMDB ID is NEVER written from memory. Every id is looked up by name
     against the sources pinned in np0 (MRCLinkdb, MetaLigand + its synonym
     lists, MetalinksDB, MEBOCOST). If the expanded name is not in any of them
     the id stays empty and the row stays flagged.

Rows where the expansion is confident get needs_review=0. Rows where the
chemistry is genuinely ambiguous keep needs_review=1 with a restated reason, so
a human still has to sign them off. Excluded rows are marked EXCLUDE.

Re-runnable: reads and rewrites the same file, and only ever touches rows whose
manual_* columns are still empty.

Usage: python3 np1b_fill_review.py
"""
from __future__ import annotations

import csv
import re
import importlib.util
import sys
from pathlib import Path

from cur_common import HERE, log

MAP = HERE / "np_composite_map.tsv"

# token -> (expanded name to look up, confidence, rationale)
#   confidence 'ok'    -> needs_review cleared
#   confidence 'check' -> id proposed, row stays flagged for sign-off
EXPANSION: dict[str, tuple[str, str, str]] = {
    "5alphaP":     ("5alpha-Dihydroprogesterone", "ok",
                    "DHRS9 reduces progesterone to 5alpha-dihydroprogesterone"),
    "Aldo":        ("Aldosterone", "ok", "CYP11B2 is aldosterone synthase"),
    "Cort":        ("Corticosterone", "ok",
                    "CYP11B1 is 11-beta-hydroxylase; resolves the CORT/cortistatin false positive"),
    "DHEAS":       ("Dehydroepiandrosterone sulfate", "ok", "SULT2B1 sulfonates DHEA"),
    "DHEAsulfate": ("Dehydroepiandrosterone sulfate", "ok", "same compound, CellPhoneDB spelling"),
    "DHT":         ("Dihydrotestosterone", "ok", "SRD5A1/2/3 are 5-alpha reductases"),
    "LTB4":        ("Leukotriene B4", "ok", "LTA4H is leukotriene A4 hydrolase"),
    "LTC4":        ("Leukotriene C4", "ok", "LTC4S is leukotriene C4 synthase"),
    "PGF2a":       ("Prostaglandin F2a", "ok", "prostaglandin F2alpha"),
    "PGI2":        ("Prostaglandin I2", "ok", "PTGIS is prostacyclin synthase"),
    "Prostacyclin": ("Prostaglandin I2", "ok", "prostacyclin is prostaglandin I2"),
    "5SHpETE":     ("5(S)-Hydroperoxyeicosatetraenoic acid", "ok",
                    "ALOX5 product; 5(S)-HPETE"),
    "atRetinoicAcid": ("all-trans-Retinoic acid", "ok", "ALDH1A1/2/3 oxidise retinal"),
    "b-PEA":       ("Phenethylamine", "ok", "beta-phenylethylamine; DDC decarboxylates phenylalanine"),

    # --- ambiguous: candidate proposed, human sign-off still required --------
    "12oxoLTB4":   ("12-Keto-leukotriene B4", "check",
                    "PTGR1 makes 12-oxo-LTB4; HMDB indexes it as 12-keto-LTB4. "
                    "Same compound on the face of it, but confirm the tautomer naming."),
    "12oxoLeukotrieneB4": ("12-Keto-leukotriene B4", "check",
                    "as 12oxoLTB4, CellPhoneDB spelling"),
    "22HC":        ("22b-Hydroxycholesterol", "check",
                    "CYP11A1 acts on 22R-hydroxycholesterol; HMDB has '22b-'. "
                    "CONFIRM 22b == 22R before accepting - stereochemistry matters here."),
    "22Hydroxycholesterol": ("22b-Hydroxycholesterol", "check",
                    "as 22HC, CellPhoneDB spelling"),
}

# Two chemicals fused into one upstream token. np2 splits on '|'.
SPLIT: dict[str, tuple[str, str]] = {
    "AldosteroneCorticosterone": ("Aldosterone", "Corticosterone"),
    "SerotoninDopamin":          ("Serotonin", "Dopamine"),
}

# Not chemicals at all - peptides and protein complexes filed under a
# non-protein annotation upstream. See PHASE_NONPEPTIDE.md defect 7.
EXCLUDE: dict[str, str] = {
    "CGA":            "peptide: thyrotropin alpha subunit (CGA+TSHB)",
    "GPHB5":          "peptide: thyrostimulin (GPHA2+GPHB5)",
    "PNOC":           "peptide: prepronociceptin",
    "b-Endorphin":    "peptide: beta-endorphin, POMC-derived",
    "hCGB1_complex":  "protein complex, not a metabolite",
    "hCGB2_complex":  "protein complex, not a metabolite",
    "hCGB3_complex":  "protein complex, not a metabolite",
    "hCGB7_complex":  "protein complex, not a metabolite",
}

# Gene symbols written wrong upstream. See defect 9.
ENZYME_FIX: dict[str, str] = {"SULT2B": "SULT2B1"}


def main() -> int:
    spec = importlib.util.spec_from_file_location("np1", HERE / "np1_parse_composites.py")
    np1 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(np1)
    idx = np1.chem_index()

    rows = list(csv.DictReader(MAP.open(encoding="utf-8"), delimiter="\t"))
    cols = list(rows[0].keys())
    for extra in ("manual_enzyme_genes", "first_pass_note"):
        if extra not in cols:
            cols.append(extra)
            for r in rows:
                r[extra] = ""

    n_ok = n_check = n_excl = n_split = n_noid = 0

    for r in rows:
        if r["needs_review"] != "1" or r["manual_chem_name"]:
            continue
        tok = r["chem_token"]

        # Enzymes the parser could not take because the symbol is wrong upstream.
        # np1 left the whole string in chem_token when the right-walk found no
        # gene, so strip the corrected enzyme off to recover the chemical too.
        for bad, good in ENZYME_FIX.items():
            if bad in r["raw_ligand"] and not r["enzyme_genes"]:
                r["manual_enzyme_genes"] = good
                r["first_pass_note"] = f"upstream symbol '{bad}' corrected to '{good}'"
                tok = re.sub(r"[_-]+(?:and_|by)?" + re.escape(bad) + r"$", "", tok)
                r["chem_token"] = tok

        if tok in EXCLUDE:
            r["manual_chem_name"] = "EXCLUDE"
            r["first_pass_note"] = EXCLUDE[tok]
            r["review_reason"] = "excluded: not a chemical ligand"
            n_excl += 1
            continue

        if tok in SPLIT:
            a, b = SPLIT[tok]
            ha, hb = idx.get(np1.norm_chem(a)), idx.get(np1.norm_chem(b))
            if ha and hb and ha[1] and hb[1]:
                r["manual_chem_name"] = f"{a}|{b}"
                r["manual_hmdb_id"] = f"{ha[1]}|{hb[1]}"
                r["needs_review"] = "0"
                r["first_pass_note"] = "fused upstream; split into two chemicals, '|' separated"
                n_split += 1
            continue

        if tok not in EXPANSION:
            continue
        name, conf, why = EXPANSION[tok]
        hit = idx.get(np1.norm_chem(name))
        r["manual_chem_name"] = name
        r["first_pass_note"] = why
        if not (hit and hit[1]):
            r["review_reason"] = f"expanded to '{name}' but no HMDB id in any pinned source"
            n_noid += 1
            continue
        r["manual_hmdb_id"] = hit[1]
        if conf == "ok":
            r["needs_review"] = "0"
            n_ok += 1
        else:
            r["review_reason"] = f"CONFIRM: {why}"
            n_check += 1

    with MAP.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    left = sum(1 for r in rows if r["needs_review"] == "1")
    log("np1b", f"resolved={n_ok} split={n_split} excluded={n_excl} "
                f"awaiting_confirmation={n_check} expanded_but_no_id={n_noid}")
    log("np1b", f"{len(rows) - left}/{len(rows)} tokens now clear, {left} still flagged")
    log("np1b", f"every HMDB id came from a np0-pinned source; none written from memory")
    return 0


if __name__ == "__main__":
    sys.exit(main())
