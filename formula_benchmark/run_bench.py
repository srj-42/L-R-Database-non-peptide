"""NPL availability-formula benchmark (2026-11-09).

Gene sets: curation.db ligand_gene_set (np23; read-only), identical for every formula.
Expression: HPA rna_tissue_consensus (nTPM) and rna_single_cell_type (nCPM), log2(x+1).
Truth: truth.py (known producer tissues / cell types, written before scoring).
Metric: per-ligand AUROC (producers vs all other samples); variants compared with the
current site formula (usage.js) by paired Wilcoxon over ligands; random-gene null.

Run:  python -I run_bench.py   (from any dir; paths are absolute)
"""
import sqlite3, zipfile, io, json, sys, collections
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import roc_auc_score

HERE = "/home/sreeja/omics/database/curation/work/benchmark/np_formula"
DB = "/home/sreeja/omics/database/curation/curation.db"
sys.path.insert(0, HERE)
from truth import TISSUE, CELLTYPE  # noqa: E402

ROLES = ["synthesis", "precursor_synthesis", "influx", "efflux", "degradation"]


def load_hpa(fname, col):
    with zipfile.ZipFile(f"{HERE}/dl/{fname}.tsv.zip") as z:
        df = pd.read_csv(z.open(f"{fname}.tsv"), sep="\t")
    m = df.pivot_table(index="Gene name", columns=df.columns[2], values=col, aggfunc="max")
    return np.log2(m.fillna(0) + 1)


def gene_sets():
    c = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    rows = c.execute("""select e.name, g.gene_symbol, g.role from ligand_gene_set g
                        join chemical_entities e using(chemical_id)
                        where g.hub_cofactor=0 and g.role!=''""").fetchall()
    gs = collections.defaultdict(lambda: collections.defaultdict(set))
    for n, g, r in rows:
        gs[n][r].add(g)
    return gs


# ---------- formulas: v = {role: matrix rows (genes x samples)} -> score per sample ----------
def agg(M, how):
    if M is None or len(M) == 0:
        return None
    if how == "mean":
        return M.mean(0)
    if how == "max":
        return M.max(0)
    if how == "gmean":
        return np.exp(np.log(M + 1e-3).mean(0))


def ours(v, how="mean", lam=1.0, mode="sub", use_eff=True, use_salv=True, comb="geo"):
    syn, pre, inf, deg, eff = (agg(v.get(r), how) for r in ROLES)
    if syn is None:
        return None
    pair = (lambda a, b: np.sqrt(a * b)) if comb == "geo" else (lambda a, b: (a + b) / 2)
    SD = pair(syn, pre) if pre is not None else syn
    SS = pair(syn, inf) if (inf is not None and use_salv) else None
    S = (SD + SS) / 2 if SS is not None else SD
    if deg is None or lam == 0:
        N = S
    elif mode == "sub":
        N = np.maximum(0, S - lam * deg)
    elif mode == "ratio":
        N = S * S / (S + deg + 1e-9)
    M = pair(N, eff) if (eff is not None and use_eff) else N
    return M


def mebocost(v):
    syn, deg = agg(v.get("synthesis"), "mean"), agg(v.get("degradation"), "mean")
    if syn is None:
        return None
    return syn - (deg if deg is not None else 0)


def cpdb(v):
    syn, eff = agg(v.get("synthesis"), "mean"), agg(v.get("efflux"), "mean")
    if syn is None:
        return None
    return np.sqrt(syn * eff) if eff is not None else syn


def mrclinkdb(v):
    syn, deg, eff = (agg(v.get(r), "mean") for r in ("synthesis", "degradation", "efflux"))
    if syn is None:
        return None
    a = syn - (deg if deg is not None else 0)
    return np.sqrt(a ** 2 + (eff if eff is not None else 0) ** 2)


def neuronchat(v):
    parts = [v[r] for r in ("synthesis", "efflux") if v.get(r) is not None and len(v[r])]
    if not parts or v.get("synthesis") is None:
        return None
    return agg(np.vstack(parts), "gmean")


VARIANTS = {
    "OURS (site formula)": lambda v: ours(v),
    "no degradation (=MetaLigand-like)": lambda v: ours(v, lam=0),
    "degradation x0.5": lambda v: ours(v, lam=0.5),
    "degradation as ratio": lambda v: ours(v, mode="ratio"),
    "no efflux gate": lambda v: ours(v, use_eff=False),
    "no salvage term": lambda v: ours(v, use_salv=False),
    "max instead of mean": lambda v: ours(v, how="max"),
    "max, no degradation": lambda v: ours(v, how="max", lam=0),
    "arithmetic (OR) combination": lambda v: ours(v, comb="arith"),
    "synthesis mean only": lambda v: agg(v.get("synthesis"), "mean"),
    "synthesis max only": lambda v: agg(v.get("synthesis"), "max"),
    "MEBOCOST (syn-deg)": mebocost,
    "CellPhoneDB v5 (sqrt syn*eff)": cpdb,
    "MRCLinkdb": mrclinkdb,
    "NeuronChat (gmean syn+eff)": neuronchat,
}


def scale(X, kind):
    if kind == "log":
        return X
    if kind == "rank":  # per-gene percentile across samples: removes gene-baseline differences
        return X.rank(axis=1, pct=True)


def run(X, truth, gs, label, rng, n_null=200):
    genes = set(X.index)
    res = collections.defaultdict(dict)
    zero_prod = collections.Counter()
    null = {}
    allgenes = np.array(sorted(genes))
    for lig, prod in truth.items():
        if lig not in gs:
            print(f"  [{label}] {lig}: no genes in ligand_gene_set, skipped")
            continue
        y = np.array([s in prod for s in X.columns])
        missing = [p for p in prod if p not in X.columns]
        if missing:
            print(f"  [{label}] {lig}: producer not in data {missing}")
        if y.sum() == 0:
            continue
        v = {r: X.loc[sorted(gs[lig][r] & genes)].values for r in ROLES if gs[lig].get(r)}
        v = {r: m for r, m in v.items() if len(m)}
        if "synthesis" not in v:
            print(f"  [{label}] {lig}: no synthesis gene measured, skipped")
            continue
        for name, f in VARIANTS.items():
            s = f(v)
            if s is None:
                continue
            res[name][lig] = roc_auc_score(y, s)
            if name.startswith("OURS") and np.all(s[y] == 0):
                zero_prod[lig] += 1
        # random-gene null for OURS: same role sizes
        sizes = {r: len(m) for r, m in v.items()}
        aucs = []
        for _ in range(n_null):
            vr = {r: X.loc[rng.choice(allgenes, k, replace=False)].values for r, k in sizes.items()}
            aucs.append(roc_auc_score(y, ours(vr)))
        null[lig] = np.array(aucs)
    return res, zero_prod, null


def summarise(res, null, label):
    base = res["OURS (site formula)"]
    ligs = sorted(base)
    out = []
    for name in VARIANTS:
        common = [l for l in ligs if l in res[name]]
        a = np.array([res[name][l] for l in common])
        b = np.array([base[l] for l in common])
        d = a - b
        p = wilcoxon(a, b).pvalue if np.any(d != 0) and len(d) > 5 else float("nan")
        out.append(dict(dataset=label, variant=name, n=len(common), mean_auc=a.mean(),
                        median_auc=np.median(a), n_auc_ge_0_8=int((a >= 0.8).sum()),
                        wins=int((d > 1e-9).sum()), losses=int((d < -1e-9).sum()),
                        p_vs_ours=p))
    nulls = np.array([np.mean(null[l]) for l in ligs])
    emp = np.array([(null[l] >= base[l]).mean() for l in ligs])
    out.append(dict(dataset=label, variant="random genes (null, OURS formula)", n=len(ligs),
                    mean_auc=nulls.mean(), median_auc=np.median(nulls), n_auc_ge_0_8=np.nan,
                    wins=np.nan, losses=np.nan, p_vs_ours=np.nan))
    return pd.DataFrame(out), dict(zip(ligs, emp))


def main():
    rng = np.random.default_rng(20261109)
    gs = gene_sets()
    tables, per_lig = [], []
    for fname, col, truth, label in [
        ("rna_tissue_consensus", "nTPM", TISSUE, "HPA tissues (50)"),
        ("rna_single_cell_type", "nCPM", CELLTYPE, "HPA cell types"),
    ]:
        X0 = load_hpa(fname, col)
        for sc in ("log", "rank"):
            X = scale(X0, sc)
            lab = f"{label} | {sc}"
            print("==", lab, X.shape)
            res, zp, null = run(X, truth, gs, lab, rng)
            t, emp = summarise(res, null, lab)
            tables.append(t)
            for l in sorted(res["OURS (site formula)"]):
                row = {"dataset": lab, "ligand": l, "producers": "; ".join(truth[l]),
                       "OURS_all_producers_zero": bool(zp.get(l)), "null_p_OURS": emp[l]}
                row.update({k: res[k].get(l) for k in VARIANTS})
                per_lig.append(row)
    T = pd.concat(tables)
    T.to_csv(f"{HERE}/summary.csv", index=False)
    pd.DataFrame(per_lig).to_csv(f"{HERE}/per_ligand.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 20, "display.max_rows", 200)
    print(T.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
