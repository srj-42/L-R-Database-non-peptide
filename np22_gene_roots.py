#!/usr/bin/env python3
"""np22_gene_roots.py - ligand gene sets rebuilt from ROOT sources, and a check of
the aggregated metabolite_enzymes rows against them (owner request 2026-11-07).

WHY. metabolite_enzymes aggregates six resources, each of which found its genes a
different way:
  MetaLigand   - Human-GEM / Mouse-GEM (Metabolic Atlas): 'produce' reactions = synthesis
                 (step 1), same for precursors (step 2), transport reactions cytosol->
                 extracellular = efflux, extracellular->cytosol for precursors = influx;
                 then manual curation (literature, TCDB). No degradation genes.
  MetalinksDB  - Recon3D + Human-GEM ('hmr') + Rhea + HMDB: producer (+1) / consumer (-1);
                 reversible reactions in both directions (transporters excepted);
                 transport direction from compartments.
  MEBOCOST     - HMDB reaction records (metabolite, reaction, gene, substrate/product).
  MRCLinkdb    - HMDB protein associations (HMDBP ids) with product/substrate flags;
                 transporters without direction.
  CellPhoneDB  - hand-curated last biosynthesis enzyme ('<ligand>_by<GENE>').
  CellChatDB   - hand-curated synthase + transporter composite per ligand.
So the roots are: genome-scale models, Rhea, HMDB (plus TCDB for transporters).
Recon3D is not used separately: Human-GEM was built by merging HMR2, Recon3D and iHsa.

WHAT THIS SCRIPT BUILDS (table ligand_gene_root, one row per gene x role x source x ref):
  human_gem    Human-GEM v2.1.1 (raw/human_gem_v2.1.1, pinned, SHA256SUMS).
               Reaction with a gene rule where the ligand is
                 product only            -> synthesis      (reversible: + degradation, flagged)
                 substrate only          -> degradation    (reversible: + synthesis, flagged)
                 on both sides, c -> e   -> efflux ; e -> c -> influx (reversible: other, flagged)
                 on both sides, other compartments -> transport_<from>_<to> (organelles)
               precursor_synthesis / precursor_influx: for each non-reversed synthesis
               reaction, its other substrates that are not currency metabolites
               (CURRENCY below); genes of reactions that make them / move them e -> c.
               The model has no secretory-vesicle compartment, so vesicle loading is
               only present where the model wrote it as c -> e (e.g. SLC18A2).
  rhea_uniprot Reviewed human UniProt CATALYTIC ACTIVITY (raw/uniprot_rhea/) joined to
               Rhea reaction sides and in/out locations (Rhea SPARQL, cached in
               work/np22/). The ligand's ChEBI id (from HMDB) is expanded to its charged
               forms (RO_0018033/RO_0018034), because Rhea uses the major species at
               pH 7.3. Direction from UniProt 'PhysiologicalDirection'; none given ->
               both roles, flagged reversible. experimental = ECO:0000269 evidence.
               Rhea in/out does not say WHICH membrane (plasma vs organelle).
  tcdb         TCDB substrate (ChEBI) x UniProt accession (work/np19 caches): role
               'transport', no direction.
  hmdb         HMDB 5.0 protein associations (work/np20/hmdb_context.json): role
               'associated', no direction (HMDB's own protein_type kept in detail).

ADDED 2026-11-08 (owner: "do what you recommend"):
  reactome     Reactome human reactions (raw/reactome_np/*_human_2026-11-08.txt): the gene's
               protein and the ligand (ChEBI incl. charged forms) are participants of the same
               curated reaction. Role unknown -> 'associated'. Reaction id in ref_id.
  go_exp       GO annotations (raw/go_np/goa_human_2026-11-08.gaf.gz, go-basic 2026-07-26)
               with experimental or curator evidence only (EXP IDA IPI IMP IGI IEP TAS IC HDA
               HMP HGI HEP HTP; NOT-qualified dropped). GO terms are matched to chemicals by
               NAME ('<X> biosynthetic process' -> synthesis, 'catabolic' -> degradation,
               'secretion'/'export' -> efflux, 'uptake'/'import' -> influx, 'transport'/
               'transporter activity' -> transport, 'metabolic process' -> associated). X must
               be the chemical's own name or a specific synonym (unique to one chemical, not a
               class word such as 'aldehyde', 'nucleotide' or 'hexose'); a synonym-matched term must have
               <= 2 sub-terms and passes nothing on. Role-specific subtypes (is_a children) of
               own-name terms inherit; 'metabolic process' terms do not. GO process annotations
               include regulators (e.g. DRD2 'arachidonate secretion'), so GO is weaker than
               Human-GEM/Rhea: it counts toward 'supported' but never sets a direction on its own.
  Support rule (owner-approved recommendation): a metabolite_enzymes row is SUPPORTED if
  Human-GEM, Rhea, Reactome or experimental GO links that gene to that chemical (or TCDB, for
  transport roles). HMDB alone does not count. Unsupported rows are kept and labelled.
  hub_cofactor: chemicals that are metabolic currency (ATP, ADP, phosphate, NAD, CoA, ...)
  or have > 300 root genes: their gene sets are not usable as an expression proxy.

AND a check table metabolite_enzyme_support: for every metabolite_enzymes row, which
roots support the same gene, and whether in the same role.

Both tables describe molecules, not ligand-receptor pairs; no verification status
depends on them. Gene symbols are normalised to current HGNC (raw/hgnc).
Usage:
  python3 np22_gene_roots.py --fetch-rhea        # SPARQL -> work/np22/*.json
  python3 np22_gene_roots.py [--commit]          # build + load (dry run without)
"""
import csv, collections, json, pickle, re, sqlite3, sys, time, urllib.parse, urllib.request
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
GEM = HERE / "raw/human_gem_v2.1.1"
GEM_VERSION = "Human-GEM v2.1.1"
UNIPROT = HERE / "raw/uniprot_rhea/uniprot_human_catalytic_2026-11-07.tsv"
WORK = HERE / "work/np22"
HMDB_CACHE = HERE / "work/np20/hmdb_context.json"
TCDB_SUB = HERE / "work/np19/tcdb_substrates.tsv"
TCDB_ACC = HERE / "work/np19/tcdb_acc2tcid.tsv"
HGNC = HERE / "raw/hgnc/hgnc_complete_set.txt"
REACT_CHEBI = HERE / "raw/reactome_np/ChEBI2Reactome_PE_Reactions_human_2026-11-08.txt"
REACT_UNI = HERE / "raw/reactome_np/UniProt2Reactome_PE_Reactions_human_2026-11-08.txt"
GO_GAF = HERE / "raw/go_np/goa_human_2026-11-08.gaf.gz"
GO_OBO = HERE / "raw/go_np/go-basic_2026-11-08.obo"
GO_EVID = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC", "HDA", "HMP", "HGI", "HEP", "HTP"}
GENERIC_SYN = {"aldehyde", "nucleotide", "dinucleotide", "acid", "alcohol", "amine", "sugar", "lipid",
               "fatty acid", "steroid", "vitamin", "ketone", "ester", "peptide", "protein", "hormone",
               "sterol", "glycoside", "nucleoside", "purine", "pyrimidine", "amino acid", "oxylipin",
               "eicosanoid", "prostaglandin", "phospholipid", "carbohydrate", "monosaccharide",
               "hexose", "d-hexose", "pentose", "aldose", "ketose", "disaccharide", "oligosaccharide",
               "polysaccharide", "bile acid", "bile salt", "catecholamine", "indoleamine", "monoamine",
               "polyamine", "nucleobase", "ribonucleotide", "deoxyribonucleotide", "glycerolipid",
               "sphingolipid", "ceramide", "neurotransmitter", "organic acid", "carboxylic acid", "cation", "anion"}
HUB_NAMES = {"Adenosine triphosphate", "ADP", "Adenosine monophosphate", "Phosphate", "Pyrophosphate",
             "NAD", "NADH", "NADP", "NADPH", "Coenzyme A", "Acetyl-CoA", "Guanosine triphosphate",
             "Guanosine diphosphate", "Uridine triphosphate", "Uridine 5'-diphosphate",
             "Cytidine triphosphate", "CDP", "Hydrogen Ion", "Water", "Oxygen", "Carbon dioxide",
             "Ammonium", "Glutathione", "S-Adenosylmethionine", "Flavin adenine dinucleotide"}
csv.field_size_limit(10**8)

CURRENCY = {"H+", "H2O", "CoA", "ATP", "Na+", "NADP+", "NADPH", "Pi", "ADP", "O2", "NAD+",
            "HCO3-", "NADH", "PPi", "AMP", "UDP", "FAD", "FADH2", "CO2", "H2O2", "chloride",
            "K+", "NH3", "GTP", "GDP", "UTP", "CTP", "GMP", "UMP", "CMP", "Fe2+", "Fe3+",
            "Ca2+", "sulfate", "PAPS", "PAP", "SAM", "SAH", "GSH", "GSSG", "acetyl-CoA",
            "L-carnitine", "O2-", "NO", "Mg2+", "Zn2+", "Cu2+", "Mn2+", "iodide", "2-oxoglutarate",
            "succinate", "fumarate", "malonyl-CoA", "tetrahydrobiopterin", "dihydrobiopterin",
            "thioredoxin", "ubiquinol", "ubiquinone", "FMN", "UDP-glucuronate", "UDP-glucose",
            "formate", "acetate", "pyruvate", "Cl-", "ferricytochrome C", "ferrocytochrome C"}

def log(*a):
    print(*a, flush=True)

# ------------------------------------------------------------------ HGNC normaliser
def hgnc_map():
    cur, alias = {}, {}
    for r in csv.DictReader(open(HGNC), delimiter="\t"):
        s = r["symbol"]
        cur[s.upper()] = s
        for col in ("prev_symbol", "alias_symbol"):
            for a in (r.get(col) or "").strip('"').split("|"):
                if a:
                    alias.setdefault(a.upper(), s)
    return lambda g: cur.get(g.upper()) or alias.get(g.upper()) or g

# ------------------------------------------------------------------ Rhea SPARQL fetch
RHEA_EP = "https://sparql.rhea-db.org/sparql?"
PFX = """PREFIX rh:<http://rdf.rhea-db.org/>
PREFIX CHEBI:<http://purl.obolibrary.org/obo/CHEBI_>
PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl:<http://www.w3.org/2002/07/owl#>
PREFIX RO:<http://purl.obolibrary.org/obo/RO_>
"""

def sparql(q):
    url = RHEA_EP + urllib.parse.urlencode({"query": PFX + q, "format": "json"})
    for attempt in range(4):
        try:
            r = urllib.request.urlopen(urllib.request.Request(
                url, headers={"Accept": "application/sparql-results+json",
                              "User-Agent": "lr-curation-script"}), timeout=180)
            return json.load(r)["results"]["bindings"]
        except Exception as e:                         # noqa: BLE001
            log("  sparql retry", attempt, e)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("Rhea SPARQL failed")

def our_chebis():
    hc = json.load(open(HMDB_CACHE))
    c = sqlite3.connect(f"file:{HERE/'curation.db'}?mode=ro", uri=True)
    out = {}
    for cid, h in c.execute("SELECT chemical_id, hmdb_id FROM chemical_entities"):
        ch = ((hc.get(h or "") or {}).get("chebi_id") or "").strip()
        if ch:
            out[cid] = ch
    return out

def fetch_rhea():
    WORK.mkdir(parents=True, exist_ok=True)
    chebis = sorted(set(our_chebis().values()))
    log("ChEBI ids to expand:", len(chebis))
    forms = collections.defaultdict(set)              # our chebi -> all forms incl itself
    for i in range(0, len(chebis), 60):
        batch = chebis[i:i + 60]
        vals = " ".join(f"CHEBI:{x}" for x in batch)
        q = f"""SELECT ?input ?ms WHERE {{ VALUES ?input {{ {vals} }}
          {{ ?ms rdfs:subClassOf ?r . ?r owl:onProperty ?p ; owl:someValuesFrom ?input .
             VALUES ?p {{ RO:0018033 RO:0018034 }} }} UNION
          {{ ?input rdfs:subClassOf ?r . ?r owl:onProperty ?p ; owl:someValuesFrom ?ms .
             VALUES ?p {{ RO:0018033 RO:0018034 }} }} }}"""
        for b in sparql(q):
            forms[b["input"]["value"].rsplit("_", 1)[1]].add(b["ms"]["value"].rsplit("_", 1)[1])
        for x in batch:
            forms[x].add(x)
    allforms = sorted({f for s in forms.values() for f in s})
    log("forms incl charged species:", len(allforms))
    parts = collections.defaultdict(list)             # chebi form -> [(master, side, loc)]
    for i in range(0, len(allforms), 40):
        vals = " ".join(f"CHEBI:{x}" for x in allforms[i:i + 40])
        q = f"""SELECT ?chebi ?rhea ?side ?loc WHERE {{ VALUES ?chebi {{ {vals} }}
          ?c rh:chebi ?chebi . ?part rh:compound ?c . ?side rh:contains ?part .
          ?rhea rh:side ?side . OPTIONAL {{ ?part rh:location ?loc }} }}"""
        for b in sparql(q):
            parts[b["chebi"]["value"].rsplit("_", 1)[1]].append(
                (b["rhea"]["value"].rsplit("/", 1)[1], b["side"]["value"][-1],
                 b.get("loc", {}).get("value", "").rsplit("/", 1)[-1]))
        log(f"  participants batch {i // 40 + 1}/{(len(allforms) + 39) // 40}")
    json.dump({k: sorted(v) for k, v in forms.items()}, open(WORK / "chebi_forms.json", "w"))
    json.dump(parts, open(WORK / "rhea_participants.json", "w"))
    log("saved", WORK / "chebi_forms.json", WORK / "rhea_participants.json")

# ------------------------------------------------------------------ Human-GEM
def load_gem():
    cache = WORK / "gem_parsed.pkl"
    if cache.exists():
        return pickle.load(open(cache, "rb"))
    y = dict(yaml.load(open(GEM / "Human-GEM.yml"), Loader=yaml.CSafeLoader))
    rx = []
    for r in y["reactions"]:
        r = dict(r)
        rx.append({"id": r["id"], "name": r.get("name", ""),
                   "mets": dict(r.get("metabolites") or []),
                   "lb": r.get("lower_bound", 0), "gpr": r.get("gene_reaction_rule") or "",
                   "subsys": "; ".join(r.get("subsystem") or [])})
    # YAML 1.1 reads the bare name NO (nitric oxide) as boolean False: turn it back
    fixb = lambda n: ("NO" if n is False else "YES" if n is True else str(n))
    mname = {dict(m)["id"][:-1]: fixb(dict(m)["name"]) for m in y["metabolites"]}
    sym = {}
    for r in csv.DictReader(open(GEM / "genes.tsv"), delimiter="\t"):
        sym[r["genes"]] = (r["geneSymbols"] or r["genes"]).split(";")[0]
    uni2sym = {}
    for r in csv.DictReader(open(GEM / "genes.tsv"), delimiter="\t"):
        for u in (r["geneUniProtID"] or "").split(";"):
            if u.strip():
                uni2sym[u.strip()] = sym[r["genes"]]
    h2m, c2m = collections.defaultdict(set), collections.defaultdict(set)
    p2m = collections.defaultdict(set)
    for r in csv.DictReader(open(GEM / "metabolites.tsv"), delimiter="\t"):
        for pc in (r["metPubChemID"] or "").split(";"):
            if pc.strip():
                p2m[pc.strip()].add(r["metsNoComp"])
        for h in (r["metHMDBID"] or "").split(";"):
            if h.strip():
                h2m[h.strip()].add(r["metsNoComp"])
        for ch in (r["metChEBIID"] or "").split(";"):
            if ch.strip():
                c2m[ch.strip().replace("CHEBI:", "")].add(r["metsNoComp"])
    out = (rx, mname, sym, uni2sym, dict(h2m), dict(c2m), dict(p2m))
    WORK.mkdir(parents=True, exist_ok=True)
    pickle.dump(out, open(cache, "wb"))
    return out

def gpr_genes(gpr):
    best = {}
    for clause in re.split(r"\bor\b", gpr.replace("(", " ").replace(")", " ")):
        gs = re.findall(r"ENSG\d+", clause)
        for g in gs:
            best[g] = best.get(g, True) and len(gs) > 1   # in_complex only if always ANDed
    return best

def gem_roles(bases, rx):
    """yield (role, reaction, reversed_flag)"""
    for r in rx:
        if not r["gpr"].strip():
            continue
        mine = {m: v for m, v in r["mets"].items() if m[:-1] in bases}
        if not mine:
            continue
        rev = r["lb"] < 0
        sub = [m for m, v in mine.items() if v < 0]
        pro = [m for m, v in mine.items() if v > 0]
        if sub and pro:
            for a in {m[-1] for m in sub}:
                for b in {m[-1] for m in pro}:
                    if a == b:
                        continue
                    role = "efflux" if (a, b) == ("c", "e") else "influx" if (a, b) == ("e", "c") \
                        else f"transport_{a}_{b}"
                    yield role, r, False
                    if rev:
                        yield ({"efflux": "influx", "influx": "efflux"}.get(role)
                               or f"transport_{b}_{a}"), r, True
        elif pro:
            yield "synthesis", r, False
            if rev:
                yield "degradation", r, True
        elif sub:
            yield "degradation", r, False
            if rev:
                yield "synthesis", r, True

# ------------------------------------------------------------------ UniProt / Rhea parse
def parse_uniprot():
    """-> list of (accession, gene, master_rhea, directions:set('LR','RL'), experimental, pmids)"""
    rows = []
    for r in csv.DictReader(open(UNIPROT), delimiter="\t"):
        acc, gene, txt = r["Entry"], (r["Gene Names (primary)"] or "").strip(), r["Catalytic activity"] or ""
        for blk in txt.split("CATALYTIC ACTIVITY:")[1:]:
            ids = [int(x) for x in re.findall(r"Rhea:RHEA:(\d+)", blk)]
            if not ids:
                continue
            master = ids[0]
            dirs = set()
            for m in re.finditer(r"PhysiologicalDirection=(left-to-right|right-to-left)", blk):
                dirs.add("LR" if m.group(1) == "left-to-right" else "RL")
            exp = "ECO:0000269" in blk
            pm = sorted(set(re.findall(r"PubMed:(\d+)", blk)))
            rows.append((acc, gene.split(";")[0], str(master), dirs, exp, pm))
    return rows

# ------------------------------------------------------------------ build
def build():
    norm = hgnc_map()
    rx, mname, sym, uni2sym, h2m, c2m, p2m = load_gem()
    lname2base = collections.defaultdict(set)
    for b_, n_ in mname.items():
        lname2base[n_.lower()].add(b_)
    name2base = collections.defaultdict(set)
    for b, n in mname.items():
        name2base[n].add(b)
    currency_bases = {b for n in CURRENCY for b in name2base.get(n, ())}
    hc = json.load(open(HMDB_CACHE))
    c = sqlite3.connect(f"file:{HERE/'curation.db'}?mode=ro", uri=True)
    chems = c.execute("SELECT chemical_id, hmdb_id, name FROM chemical_entities").fetchall()
    cid_pc = dict(c.execute("SELECT chemical_id, pubchem_cid FROM chemical_entities"))
    map_how = collections.Counter()
    rows = []   # (chemical_id, gene, role, source, ref_id, detail, reversible, in_complex, experimental, version)
    mapped = 0

    # index reactions by metabolite base for precursor lookups
    produced_by = collections.defaultdict(list)
    imported_by = collections.defaultdict(list)
    for r in rx:
        if not r["gpr"].strip():
            continue
        subs = {m[:-1] for m, v in r["mets"].items() if v < 0}
        for m, v in r["mets"].items():
            if v > 0 and m[:-1] not in subs:
                produced_by[m[:-1]].append(r)
            if v > 0 and m.endswith("c") and (m[:-1] + "e") in r["mets"] and r["mets"][m[:-1] + "e"] < 0:
                imported_by[m[:-1]].append(r)

    def add_gem(cid, role, r, rev, detail_extra=""):
        for g, cx in gpr_genes(r["gpr"]).items():
            rows.append((cid, norm(sym.get(g, g)), role, "human_gem", r["id"],
                         (r["name"] + (" | " + r["subsys"] if r["subsys"] else "") + detail_extra)[:300],
                         int(rev), int(cx), 0, GEM_VERSION))

    for cid, h, nm in chems:
        bases = set(h2m.get(h or "", set()))
        ch = ((hc.get(h or "") or {}).get("chebi_id") or "")
        how = "hmdb" if bases else ""
        if not bases and ch:
            bases = set(c2m.get(ch, set())); how = "chebi" if bases else ""
        if not bases and cid_pc.get(cid):
            bases = set(p2m.get(str(cid_pc[cid]).split(".")[0], set())); how = "pubchem" if bases else ""
        if not bases and nm:
            bases = set(lname2base.get(nm.lower(), set())); how = "name" if bases else ""
        if not bases:
            continue
        map_how[how] += 1
        mapped += 1
        precursors = set()
        for role, r, rev in gem_roles(bases, rx):
            add_gem(cid, role, r, rev)
            if role == "synthesis" and not rev:
                for m, v in r["mets"].items():
                    b = m[:-1]
                    if v < 0 and b not in bases and b not in currency_bases:
                        precursors.add(b)
        for p in precursors:
            for r in produced_by.get(p, []):
                if r["lb"] < 0:
                    continue
                add_gem(cid, "precursor_synthesis", r, False, f" | precursor: {mname.get(p, p)}")
            for r in imported_by.get(p, []):
                add_gem(cid, "precursor_influx", r, r["lb"] < 0, f" | precursor: {mname.get(p, p)}")
    log("human_gem: chemicals mapped to the model", mapped, "of", len(chems), dict(map_how))

    # ---- Rhea / UniProt
    forms = json.load(open(WORK / "chebi_forms.json"))
    parts = json.load(open(WORK / "rhea_participants.json"))
    sides_of = collections.defaultdict(lambda: collections.defaultdict(set))   # chem -> master -> {(side,loc)}
    chem_chebi = {}
    for cid, h, nm in chems:
        ch = ((hc.get(h or "") or {}).get("chebi_id") or "")
        if ch:
            chem_chebi[cid] = ch
            for f in forms.get(ch, [ch]):
                for master, side, loc in parts.get(f, []):
                    sides_of[cid][master].add((side, loc))
    up = parse_uniprot()
    by_master = collections.defaultdict(list)
    for acc, gene, master, dirs, exp, pm in up:
        by_master[master].append((acc, gene, dirs, exp, pm))
    n_rhea = 0
    for cid, masters in sides_of.items():
        for master, sl in masters.items():
            for acc, gene, dirs, exp, pm in by_master.get(master, []):
                if not gene:
                    continue
                sides = {s for s, _ in sl}
                d = ("LR" if dirs == {"LR"} else "RL" if dirs == {"RL"} else None)
                undirected = d is None
                if sides == {"L", "R"}:
                    locL = {l for s, l in sl if s == "L"}
                    locR = {l for s, l in sl if s == "R"}
                    if "Out" in locL and "In" in locR:
                        fwd = "influx"
                    elif "In" in locL and "Out" in locR:
                        fwd = "efflux"
                    else:
                        continue
                    back = {"influx": "efflux", "efflux": "influx"}[fwd]
                    roles = [(fwd, 0)] if d == "LR" else [(back, 0)] if d == "RL" else [(fwd, 1), (back, 1)]
                else:
                    right = sides == {"R"}
                    fwd = "synthesis" if right else "degradation"
                    back = "degradation" if right else "synthesis"
                    roles = [(fwd, 0)] if d == "LR" else [(back, 0)] if d == "RL" else [(fwd, 1), (back, 1)]
                for role, flag in roles:
                    rows.append((cid, norm(gene), role, "rhea_uniprot", f"RHEA:{master}",
                                 f"UniProt {acc}" + ("; PMIDs " + ",".join(pm[:5]) if pm else "")
                                 + ("; direction not stated" if undirected else ""),
                                 flag, 0, int(exp), "UniProt 2026-11-07 / Rhea SPARQL 2026-11-07"))
                    n_rhea += 1
    log("rhea_uniprot rows:", n_rhea)

    # ---- TCDB
    acc2tc = collections.defaultdict(set)
    for line in open(TCDB_ACC):
        a, t = line.rstrip("\n").split("\t")[:2]
        acc2tc[t].add(a)
    uni2gene = dict(uni2sym)
    for v in hc.values():
        for g, u, ptype, hp in v.get("proteins", []):
            if u and g:
                uni2gene.setdefault(u, g)
    for acc, gene, *_ in up:
        if gene:
            uni2gene.setdefault(acc, gene)
    n_tc = 0
    sub2tc = collections.defaultdict(set)
    for line in open(TCDB_SUB):
        tc, subs = (line.rstrip("\n").split("\t") + [""])[:2]
        for s in subs.split("|"):
            m = re.match(r"CHEBI:(\d+)", s)
            if m:
                sub2tc[m.group(1)].add(tc)
    for cid, ch in chem_chebi.items():
        for f in forms.get(ch, [ch]):
            for tc in sub2tc.get(f, ()):
                for a in acc2tc.get(tc, ()):
                    g = uni2gene.get(a)
                    if g:
                        rows.append((cid, norm(g), "transport", "tcdb", f"TC:{tc}", f"UniProt {a}",
                                     1, 0, 0, "TCDB (np19 cache 2026-10-31)"))
                        n_tc += 1
    log("tcdb rows:", n_tc)

    # ---- HMDB protein associations
    n_h = 0
    for cid, h, nm in chems:
        for g, u, ptype, hp in (hc.get(h or "") or {}).get("proteins", []):
            if g:
                rows.append((cid, norm(g), "associated", "hmdb", hp, f"{ptype or 'Unknown'}; UniProt {u}",
                             1, 0, 0, "HMDB 5.0"))
                n_h += 1
    log("hmdb rows:", n_h)
    return rows

def gene_uniprot_map(c, hc):
    g2u = collections.defaultdict(set)
    for r in csv.DictReader(open(GEM / "genes.tsv"), delimiter="\t"):
        sy = (r["geneSymbols"] or "").split(";")[0]
        for u in (r["geneUniProtID"] or "").split(";"):
            if u:
                g2u[sy].add(u)
    for v in hc.values():
        for g, u, t, hp in v.get("proteins", []):
            if g and u:
                g2u[g].add(u)
    for r in csv.DictReader(open(UNIPROT), delimiter="\t"):
        g = (r["Gene Names (primary)"] or "").split(";")[0]
        if g:
            g2u[g].add(r["Entry"])
    for sy, u in c.execute("SELECT symbol, uniprot FROM entities WHERE uniprot IS NOT NULL AND uniprot!=''"):
        g2u[sy].add(u)
    return g2u

def build_reactome(c, hc, forms, norm):
    """chemical x gene sharing a curated human Reactome reaction -> role 'associated'."""
    ch2rx = collections.defaultdict(set)
    rxname = {}
    for line in open(REACT_CHEBI):
        p = line.rstrip("\n").split("\t")
        ch2rx[p[0]].add(p[3]); rxname[p[3]] = p[5]
    u2rx = collections.defaultdict(set)
    for line in open(REACT_UNI):
        p = line.rstrip("\n").split("\t")
        u2rx[p[0].split("-")[0]].add(p[3]); rxname.setdefault(p[3], p[5])
    rx2u = collections.defaultdict(set)
    for u, rs in u2rx.items():
        for r in rs:
            rx2u[r].add(u)
    u2g = {}
    for g, us in gene_uniprot_map(c, hc).items():
        for u in us:
            u2g.setdefault(u, g)
    rows = []
    for cid, h in c.execute("SELECT chemical_id, hmdb_id FROM chemical_entities"):
        ch = ((hc.get(h or "") or {}).get("chebi_id") or "")
        if not ch:
            continue
        rxs = set()
        for f in forms.get(ch, [ch]):
            rxs |= ch2rx.get(f, set())
        for r in rxs:
            for u in rx2u.get(r, ()):
                g = u2g.get(u)
                if g:
                    rows.append((cid, norm(g), "associated", "reactome", r, rxname.get(r, "")[:250],
                                 1, 0, 0, "Reactome download 2026-11-08 (human)"))
    log("reactome rows:", len(rows))
    return rows

GO_PAT = [(r"^(.+) biosynthetic process$", "synthesis"), (r"^(.+) catabolic process$", "degradation"),
          (r"^(.+) secretion$", "efflux"), (r"^(.+) export(?: across plasma membrane| from cell)?$", "efflux"),
          (r"^(.+) uptake$", "influx"), (r"^(.+) import(?: across plasma membrane| into cell)?$", "influx"),
          (r"^(.+) transmembrane transport$", "transport"), (r"^(.+) transport$", "transport"),
          (r"^(.+) transmembrane transporter activity$", "transport"), (r"^(.+) metabolic process$", "associated")]

def build_go(c, norm):
    import gzip
    terms, parents, cur = {}, collections.defaultdict(set), None
    for line in open(GO_OBO):
        line = line.rstrip("\n")
        if line == "[Term]":
            cur = {}; continue
        if line.startswith("["):
            cur = None; continue
        if cur is None:
            continue
        if line.startswith("id: "):
            cur["id"] = line[4:]; terms[cur["id"]] = cur
        elif line.startswith("name: "):
            cur["name"] = line[6:]
        elif line.startswith("is_obsolete: true"):
            cur["obs"] = True
        elif line.startswith("is_a: "):
            parents[cur["id"]].add(line[6:16])
    # chemical name keys: own name always; synonyms only if specific
    own, syn_map = collections.defaultdict(set), collections.defaultdict(set)
    for cid, nm, syn in c.execute("SELECT chemical_id, name, synonyms FROM chemical_entities"):
        if nm:
            own[nm.lower().strip()].add(cid)
        for x in (json.loads(syn) if syn else []):
            x = (x or "").lower().strip()
            if len(x) >= 5 and not x.isdigit() and x not in GENERIC_SYN and re.search(r"[a-z]{3}", x):
                syn_map[x].add(cid)
    keys = dict(own)
    via_syn = set()
    for x, cids in syn_map.items():
        if len(cids) == 1 and x not in keys:
            keys[x] = cids
            via_syn.add(x)
    for k in list(keys):
        if k.startswith("l-") and k[2:] not in keys:
            keys[k[2:]] = keys[k]
    matched = {}
    for t in terms.values():
        if t.get("obs") or "name" not in t:
            continue
        for pat, role in GO_PAT:
            m = re.match(pat, t["name"])
            if m and m.group(1).lower() in keys:
                matched[t["id"]] = (keys[m.group(1).lower()], role, t["name"], m.group(1).lower() in via_syn); break
    kids = collections.defaultdict(set)
    for ch, ps in parents.items():
        for p in ps:
            kids[p].add(ch)
    # a term matched only through a synonym must look chemical-specific: no more than two sub-terms
    matched = {k: v for k, v in matched.items() if not v[3] or len(kids.get(k, ())) <= 2}
    full = collections.defaultdict(list)
    for tid, (cids, role, nm, by_syn) in matched.items():
        full[tid].append((cids, role, nm))
        if role == "associated" or by_syn:
            continue
        stack, seen = list(kids.get(tid, ())), set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            full[x].append((cids, role, nm))
            stack += list(kids.get(x, ()))
    rows = []
    for line in gzip.open(GO_GAF, "rt"):
        if line[0] == "!":
            continue
        p = line.split("\t")
        if p[6] not in GO_EVID or "NOT" in p[3]:
            continue
        for cids, role, nm in full.get(p[4], ()):
            for cid in cids:
                rows.append((cid, norm(p[2]), role, "go_exp", p[4],
                             f"{terms.get(p[4], {}).get('name', '')} ({p[6]}; {p[5]})"[:250],
                             int(role in ("transport", "associated")), 0, 1, "GO 2026-07-26 / GOA human 2026-11-08"))
    log("go_exp rows:", len(rows), "| GO terms matched to chemicals:", len(matched))
    return rows

def hub_chemicals(c, root_rows):
    n = collections.Counter()
    seen = set()
    for cid, g, role, src, *_ in root_rows:
        if src in ("human_gem", "rhea_uniprot") and role in ("synthesis", "degradation") and (cid, g) not in seen:
            seen.add((cid, g)); n[cid] += 1
    hubs = {cid for cid, k in n.items() if k > 300}
    hubs |= {cid for (cid,) in c.execute(
        f"SELECT chemical_id FROM chemical_entities WHERE name IN ({','.join('?' * len(HUB_NAMES))})",
        sorted(HUB_NAMES))}
    return hubs

def support(c, root_rows):
    """metabolite_enzymes row -> which sources back it, and supported / unsupported."""
    norm = hgnc_map()
    idx = collections.defaultdict(lambda: collections.defaultdict(set))   # (cid,gene) -> source -> roles
    for cid, g, role, src, *_ in root_rows:
        idx[(cid, g)][src].add(role)
    hubs = hub_chemicals(c, root_rows)
    out = []
    for me_id, cid, gene, role in c.execute("SELECT id, chemical_id, gene_symbol, role FROM metabolite_enzymes"):
        g = norm(gene or "")
        s = idx.get((cid, g), {})
        gem, rh, go = s.get("human_gem", set()), s.get("rhea_uniprot", set()), s.get("go_exp", set())
        same = lambda roles: role in roles or (role == "influx" and "precursor_influx" in roles) \
            or (role in ("efflux", "influx") and "transport" in roles)
        gem_same, rh_same, go_same = same(gem), same(rh), same(go)
        tc = int(bool(s.get("tcdb")) and role in ("efflux", "influx"))
        rc = int(bool(s.get("reactome")))
        hm = int(bool(s.get("hmdb")))
        if gem_same or rh_same or go_same or tc:
            status = "root_same_role"
        elif gem or rh or go:
            status = "root_other_role"
        elif rc:
            status = "reactome_only"
        elif hm:
            status = "hmdb_only"
        else:
            status = "no_root"
        supported = int(status in ("root_same_role", "root_other_role", "reactome_only"))
        out.append((me_id, cid, g, role, int(gem_same), int(bool(gem)), int(rh_same), int(bool(rh)), tc, hm,
                    rc, int(go_same), int(bool(go)), supported, int(cid in hubs), status))
    return out

def main():
    if "--fetch-rhea" in sys.argv:
        fetch_rhea()
        return
    rows = build()
    c = sqlite3.connect(HERE / "curation.db")
    hc = json.load(open(HMDB_CACHE))
    forms = json.load(open(WORK / "chebi_forms.json"))
    norm = hgnc_map()
    rows += build_reactome(c, hc, forms, norm)
    rows += build_go(c, norm)
    rows = list(dict.fromkeys(rows))
    c.execute("""CREATE TABLE IF NOT EXISTS ligand_gene_root(
        chemical_id TEXT NOT NULL, gene_symbol TEXT NOT NULL, role TEXT NOT NULL, source TEXT NOT NULL,
        ref_id TEXT NOT NULL, detail TEXT, reversible INTEGER NOT NULL DEFAULT 0,
        in_complex INTEGER NOT NULL DEFAULT 0, experimental INTEGER NOT NULL DEFAULT 0,
        source_version TEXT, PRIMARY KEY (chemical_id, gene_symbol, role, source, ref_id))""")
    # metabolite_enzyme_support is this script's own derived table: rebuilt with the 2026-11-08 columns
    c.execute("DROP TABLE IF EXISTS metabolite_enzyme_support")
    c.execute("""CREATE TABLE metabolite_enzyme_support(
        me_id INTEGER PRIMARY KEY, chemical_id TEXT, gene_symbol_hgnc TEXT, role TEXT,
        gem_same_role INTEGER, gem_any_role INTEGER, rhea_same_role INTEGER, rhea_any_role INTEGER,
        tcdb_transport INTEGER, hmdb_associated INTEGER, reactome_shared_reaction INTEGER,
        go_same_role INTEGER, go_any_role INTEGER, supported INTEGER, hub_cofactor INTEGER, status TEXT)""")
    c.execute("DELETE FROM ligand_gene_root")
    c.executemany("INSERT OR IGNORE INTO ligand_gene_root VALUES (?,?,?,?,?,?,?,?,?,?)", rows)
    sup = support(c, rows)
    c.executemany("INSERT INTO metabolite_enzyme_support VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", sup)
    q = lambda s: c.execute(s).fetchall()
    log("\nligand_gene_root by source x role:")
    for r in q("SELECT source, role, COUNT(*), COUNT(DISTINCT chemical_id), COUNT(DISTINCT gene_symbol) "
               "FROM ligand_gene_root GROUP BY 1,2 ORDER BY 1,3 DESC"):
        log("  ", r)
    log("chemicals with any root gene:", q("SELECT COUNT(DISTINCT chemical_id) FROM ligand_gene_root")[0][0])
    log("\nmetabolite_enzymes support status (all rows):", q("SELECT status, COUNT(*) FROM metabolite_enzyme_support GROUP BY 1 ORDER BY 2 DESC"))
    log("supported / unsupported:", q("SELECT supported, COUNT(*) FROM metabolite_enzyme_support GROUP BY 1"))
    log("hub-cofactor chemicals:", q("SELECT COUNT(DISTINCT chemical_id), COUNT(*) FROM metabolite_enzyme_support WHERE hub_cofactor=1"))
    log("by source:")
    for r in q("""SELECT m.source_id, s.status, COUNT(*) FROM metabolite_enzyme_support s
                  JOIN metabolite_enzymes m ON m.id = s.me_id GROUP BY 1,2 ORDER BY 1,3 DESC"""):
        log("  ", r)
    if "--commit" in sys.argv:
        c.commit()
        log("COMMITTED")
    else:
        c.rollback()
        log("dry run - nothing written (use --commit)")

if __name__ == "__main__":
    main()
