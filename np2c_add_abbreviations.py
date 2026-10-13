#!/usr/bin/env python3
"""Add standard, unambiguous abbreviations to chemical_entities.synonyms.

Why this exists. np9_verify.py's chemical gate requires the name a paper uses to be a
RECORDED name or synonym of that chemical, which is the right rule: it stops a quote
about one molecule being filed under another. But the synonym lists in this build came
from upstream name fields and are both incomplete and partly corrupted (a name was
comma-split, so 'Leukotriene D4' carries fragments like '11Z' and '14-tetraenoic acid'
while lacking 'LTD4'). 'Leukotriene C4' carries no synonyms at all.

The lipid-mediator literature almost never spells these out: abstracts write LTC(4),
LTD(4), PGE(2). So correct confirmations were being refused on vocabulary, and the
tempting workaround - loosening the gate, or filing the quote under the duplicate
name-keyed entity that happens to be called 'LTC4' - would damage either the gate or
chemical identity. Completing the synonym list is the honest fix, and it is what the
gate's own refusal message tells the curator to do.

Only universal abbreviations go in, each tied to one molecule by its HMDB id, which is
asserted below before anything is written. Nothing is removed; the corrupted fragments
are left alone because other joins were built against them.

Usage: python3 np2c_add_abbreviations.py [--commit]
"""
from __future__ import annotations

import json
import re
import sys

from cur_common import connect, log, now

# hmdb_id -> (expected name fragment, abbreviations to add)
ADD = {
    "HMDB0001198": ("leukotriene c4", ["LTC4", "LTC(4)"]),
    "HMDB0003080": ("leukotriene d4", ["LTD4", "LTD(4)"]),
    "HMDB0002200": ("leukotriene e4", ["LTE4", "LTE(4)"]),
    "HMDB0001085": ("leukotriene b4", ["LTB4", "LTB(4)"]),
    # ids verified against chemical_entities, not recalled - a first draft had
    # HMDB0001335 as prostaglandin E2 when it is prostaglandin I2, and the assertion
    # below refused the write rather than mislabelling the molecule.
    "HMDB0001220": ("prostaglandin e2", ["PGE2", "PGE(2)"]),
    "HMDB0001403": ("prostaglandin d2", ["PGD2", "PGD(2)"]),
    "HMDB0001139": ("prostaglandin f2", ["PGF2alpha", "PGF(2alpha)", "PGF2a"]),
    "HMDB0001335": ("prostaglandin i2", ["PGI2", "prostacyclin"]),
    "HMDB0001452": ("thromboxane a2", ["TXA2", "TxA2"]),
    "HMDB0000277": ("sphingosine 1-phosphate", ["S1P"]),
}

# Abbreviations that cannot be derived from the name. Each is the form the primary
# literature actually uses, keyed by the chemical's own name so the match is exact.
BY_NAME = {
    "3-Hydroxybutyric acid": ["\u03b2-hydroxybutyrate"],
    "phosphatidyl (3,5) inositol bisphosphate": ["PI(3,5)P2", "PtdIns(3,5)P2", "phosphatidylinositol 3,5-bisphosphate"],
    "L-Lactic acid": ["lactate", "L-lactate", "lactic acid"],
    "3-Hydroxyoctanoic acid": ["3-OH-octanoic acid", "3-hydroxyoctanoate"],
    "Geranylgeranyl-PP": ["geranylgeranyl pyrophosphate", "geranylgeranyl-PP"],
    "b-D-hydroxybutyric acid": ["(D)-beta-Hydroxybutyrate", "D-beta-hydroxybutyrate", "D-3-hydroxybutyrate"],
    "Adenosine diphosphate ribose": ["ADPR", "ADP-ribose", "ADP ribose"],
    "Phosphoserine": ["L-serine-O-phosphate", "O-phospho-L-serine", "L-SOP"],
    "15-Deoxy-d-12,14-PGJ2": ["15dPGJ2", "15-deoxy-delta12,14-prostaglandin J2", "15-deoxy-delta 12,14-PGJ2", "15-deoxy-delta12,14-PGJ2"],
    "(R)-5-Diphosphomevalonic acid": ["mevalonate 5-diphosphate", "mevalonate-5-diphosphate", "(R)-5-diphosphomevalonate"],
    "lysophosphatidylserine": ["lysoPS", "lysophosphatidyl-L-serine"],
    "9-cis-Retinoic acid": ["9-cis-RA"],
    "all-trans-Retinoic acid": ["t-RA", "all-trans-RA"],
    "LysoPA(18:1(9Z)/0:0)": ["LPA (oleoyl)", "oleoyl-LPA", "18:1 LPA", "oleoyl-L-alpha-lysophosphatidic acid", "1-oleoyl-LPA"],
    "12S-HHT": ["12-HHT", "12-hydroxyheptadecatrienoic acid", "12(S)-hydroxyheptadeca-5Z, 8E, 10E-trienoic acid"],
    "Adenosine monophosphate": ["AMP"],
    "Sphinganine": ["dihydro-D-erythro-sphingosine", "dihydrosphingosine"],
    "Oxoglutaric acid": ["alpha-ketoglutarate", "alpha-ketoglutaric acid", "2-oxoglutarate"],
    "Sphinganine 1-phosphate": ["dihydrosphingosine 1-phosphate", "dihydro-sphingosine 1-phosphate", "dihydro-S1P"],
    "Prostaglandin F2a": ["PGF2"],
    "Prostaglandin J2": ["PGJ2", "PGJ(2)"],
    "Prostaglandin E1": ["PGE1", "PGE(1)"],
    "Prostaglandin F1a": ["PGF1a", "PGF1alpha"],
    "15S-HETE": ["15(S)-hydroxyeicosatetraenoic acid", "15(S)-HETE"],
    "12-HETE": ["12(S)-hydroxyeicosatetraenoic acid", "12(S)-HETE"],
    "Calcium": ["Ca2+", "Ca(2+)"],
    "Magnesium": ["Mg2+", "Mg(2+)"],
    "Sodium": ["Na+", "Na(+)"],
    "Zinc": ["Zn2+", "Zn(2+)"],
    "LPA": ["lysophosphatidic acid"],
    "Hydrogen Ion": ["proton", "protons"],
    "Isopentenyl pyrophosphate": ["IPP", "isopentenyl-diphosphate", "isopentenyl diphosphate"],
    "Dimethylallylpyrophosphate": ["DMAPP", "dimethylallyl-diphosphate",
                                   "3,3-dimethylallyl-diphosphate"],
    "Geranyl-PP": ["GPP", "geranyl-diphosphate", "geranyl diphosphate"],
    "Geranylgeranyl-PP": ["GGPP", "geranylgeranyl-diphosphate", "geranylgeranyl diphosphate"],
    "Farnesyl pyrophosphate": ["FPP", "farnesyl-diphosphate", "farnesyl diphosphate"],
}

# Standard three-letter codes. Added for the L- form only: a paper writing 'Lys' in an
# amino-acid series means the L-enantiomer, and attaching 'Lys' to a D- entity would be
# the stereoisomer error this build already has too much of.
AA3 = {
    "L-Alanine": "Ala", "L-Arginine": "Arg", "L-Asparagine": "Asn", "L-Aspartic acid": "Asp",
    "L-Cysteine": "Cys", "L-Glutamic acid": "Glu", "L-Glutamine": "Gln", "Glycine": "Gly",
    "L-Histidine": "His", "L-Isoleucine": "Ile", "L-Leucine": "Leu", "L-Lysine": "Lys",
    "L-Methionine": "Met", "L-Phenylalanine": "Phe", "L-Proline": "Pro", "L-Serine": "Ser",
    "L-Threonine": "Thr", "L-Tryptophan": "Trp", "L-Tyrosine": "Tyr", "L-Valine": "Val",
}


def derived(name: str) -> list[str]:
    """Only naming-convention swaps that cannot produce a non-word.

    A first version also turned any '...ic acid' into '...ate'. It produced real forms
    (Propionic acid -> Propionate) but also 'propanate', '3-Mercaptopropanate' and
    'D-Tartarate', because English carboxylate naming is irregular. Junk synonyms are not
    harmless here: they are what the chemical-identity joins run on, and this build
    already carries comma-split fragments doing damage. So the generated set is limited to
    the pyrophosphate/diphosphate swap, which is purely conventional, and everything else
    is hand-verified in BY_NAME / CARBOXYLATE / AA3 as the audit actually needs it.
    """
    out, n = [], (name or "").strip()
    if "pyrophosphate" in n.lower():
        out.append(re.sub("pyrophosphate", "diphosphate", n, flags=re.I))
    if n in AA3:
        out += [AA3[n], "L-" + AA3[n]]
        # Papers almost never write the L- prefix: "glutamine", "leucine", "serine".
        # norm("L-Glutamine") = "lglutamine" is not a substring of "glutamine", so the
        # name gate rejected correct confirmations. Safe because L- is dropped only for
        # the proteinogenic amino acids listed in AA3, where "serine" in a paper means
        # the L-enantiomer unless it says D-.
        if n.startswith("L-"):
            out.append(n[2:])
    return list(dict.fromkeys(out))


# Carboxylate forms the literature uses, verified one by one against the acid they belong
# to. Short-chain fatty acid papers write the salt form almost exclusively.
CARBOXYLATE = {
    "L-Glutamic acid": ["glutamate", "L-glutamate"],
    "L-Aspartic acid": ["aspartate", "L-aspartate"],
    "Acetic acid": ["acetate"],
    "Propionic acid": ["propionate"],
    "propanoic acid": ["propionate", "propanoate"],
    "Butyric acid": ["butyrate"],
    "Isobutyric acid": ["isobutyrate"],
    "Valeric acid": ["valerate", "pentanoate"],
    "Formic acid": ["formate"],
    "Lactic acid": ["lactate"],
    "Citric acid": ["citrate"],
    "Succinic acid": ["succinate"],
    "Fumaric acid": ["fumarate"],
    "Malic acid": ["malate"],
    "Pyruvic acid": ["pyruvate"],
    "alpha-Ketoglutaric acid": ["alpha-ketoglutarate", "2-oxoglutarate"],
    "3-Hydroxybutyric acid": ["3-hydroxybutyrate", "beta-hydroxybutyrate"],
    "Acetoacetic acid": ["acetoacetate"],
}


def main() -> int:
    commit = "--commit" in sys.argv
    c = connect()
    added = skipped = 0
    # Which normalised names already belong to some chemical - a synonym must never be
    # attached to a second molecule, or a later name-based join silently merges two.
    owner: dict[str, set[str]] = {}
    for r in c.execute("SELECT chemical_id, name, synonyms FROM chemical_entities"):
        try:
            syn = json.loads(r["synonyms"] or "[]")
        except Exception:                                    # noqa: BLE001
            syn = []
        for x in [r["name"]] + [s for s in syn if isinstance(s, str)]:
            k = re.sub(r"[^a-z0-9]", "", (x or "").lower())
            if k:
                owner.setdefault(k, set()).add(r["chemical_id"])

    plan: dict[str, list[str]] = {}
    for r in c.execute("SELECT chemical_id, name FROM chemical_entities"):
        want = (list(BY_NAME.get(r["name"], []))
                + list(CARBOXYLATE.get(r["name"], []))
                + derived(r["name"]))
        if want:
            plan[r["chemical_id"]] = want
    for cid, want in plan.items():
        r = c.execute("SELECT chemical_id, name, synonyms FROM chemical_entities "
                      "WHERE chemical_id=?", (cid,)).fetchone()
        try:
            syn = json.loads(r["synonyms"] or "[]")
            syn = syn if isinstance(syn, list) else []
        except Exception:                                    # noqa: BLE001
            syn = []
        have = {s.lower() for s in syn if isinstance(s, str)}
        new = []
        for a in want:
            if a.lower() in have:
                continue
            k = re.sub(r"[^a-z0-9]", "", a.lower())
            held = owner.get(k, set()) - {cid}
            if held:
                log("np2c", f"  SKIP {a!r} for {r['name']!r}: already a name of {sorted(held)}")
                continue
            new.append(a)
            owner.setdefault(k, set()).add(cid)
        if new:
            c.execute("UPDATE chemical_entities SET synonyms=? WHERE chemical_id=?",
                      (json.dumps(syn + new), cid))
            log("np2c", f"  {cid} {r['name'][:30]:30s} += {new}")
            added += len(new)

    for hmdb, (expect, abbrevs) in ADD.items():
        r = c.execute("SELECT chemical_id, name, synonyms FROM chemical_entities "
                      "WHERE hmdb_id = ?", (hmdb,)).fetchone()
        if r is None:
            log("np2c", f"  {hmdb}: no such chemical, skipped")
            skipped += 1
            continue
        # assert the id really is the molecule these abbreviations belong to
        if expect not in (r["name"] or "").lower():
            log("np2c", f"  ABORT: {hmdb} is {r['name']!r}, expected a name containing "
                        f"{expect!r} - refusing to attach {abbrevs}")
            return 1
        try:
            syn = json.loads(r["synonyms"] or "[]")
            syn = syn if isinstance(syn, list) else []
        except Exception:                                    # noqa: BLE001
            syn = []
        have = {s.lower() for s in syn if isinstance(s, str)}
        new = [a for a in abbrevs if a.lower() not in have]
        if not new:
            skipped += 1
            continue
        c.execute("UPDATE chemical_entities SET synonyms = ? WHERE chemical_id = ?",
                  (json.dumps(syn + new), r["chemical_id"]))
        log("np2c", f"  {r['chemical_id']} {r['name'][:34]:34s} += {new}")
        added += len(new)
    log("np2c", f"added {added} abbreviation(s); {skipped} chemical(s) already complete")
    (c.commit if commit else c.rollback)()
    log("np2c", "COMMITTED" if commit else "DRY RUN - rerun with --commit to write")
    return 0


if __name__ == "__main__":
    sys.exit(main())
