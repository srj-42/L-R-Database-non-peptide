-- Phase NP schema. Additive only: no existing table is altered except
-- `sources`, which gains `shippable` in np0_acquire.py.
--
-- Chemicals deliberately do NOT go into `entities`. A metabolite has no UniProt
-- accession and no gene symbol, and forcing it into a gene-symbol column is the
-- defect PROJECT_REPORT.md section 7 objects to. Receptors keep using `entities`.

CREATE TABLE IF NOT EXISTS chemical_entities (
    chemical_id      TEXT PRIMARY KEY,          -- internal, stable: CHEM00001
    name_key         TEXT NOT NULL UNIQUE,      -- normalised name; dedup key when hmdb_id is NULL
    hmdb_id          TEXT UNIQUE,               -- PRIMARY external key (decision 4); NULL allowed
    chebi_id         TEXT,                      -- deliberately NULL until a mapping pass
    pubchem_cid      TEXT,
    pubchem_sid      TEXT,
    name             TEXT NOT NULL,
    synonyms         TEXT,                      -- JSON list
    smiles           TEXT,
    formula          TEXT,
    -- HMDB taxonomy, from MRCLinkdb. Kept separate from the ChEBI ontology
    -- below: they are different vocabularies and merging them manufactured 213
    -- spurious conflicts on the first run.
    kingdom          TEXT,
    super_class      TEXT,
    class            TEXT,
    -- ChEBI ontology, from MetaLigand's Ligand_list (Class_source = ChEBI).
    chebi_super_class TEXT,
    chebi_class      TEXT,
    metabolite_subclass TEXT,
    molecular_class  TEXT,                      -- PROVISIONAL, see np2 derive_class()
    is_neurotransmitter INTEGER DEFAULT 0,
    n_sources        INTEGER NOT NULL DEFAULT 0,
    sources          TEXT NOT NULL,             -- JSON list of source_id
    shippable        INTEGER NOT NULL DEFAULT 1 -- 0 when only a non-shippable source has it
);

CREATE INDEX IF NOT EXISTS ix_chem_hmdb ON chemical_entities(hmdb_id);

CREATE TABLE IF NOT EXISTS chemical_interactions (
    chem_interaction_id TEXT PRIMARY KEY,
    chemical_id         TEXT NOT NULL REFERENCES chemical_entities,
    receptor_entity_id  TEXT NOT NULL REFERENCES entities,
    species             TEXT NOT NULL,
    complex_id          TEXT,
    complex_role        TEXT,
    receptor_type       TEXT,
    n_sources           INTEGER NOT NULL DEFAULT 0,
    sources             TEXT NOT NULL,
    has_citation        INTEGER NOT NULL DEFAULT 0,
    n_pmids             INTEGER NOT NULL DEFAULT 0,
    verification_status TEXT DEFAULT 'unverified',
    shippable           INTEGER NOT NULL DEFAULT 1,
    -- NO confidence_score column in v1 (decision 9)
    UNIQUE (chemical_id, receptor_entity_id, species)
);

CREATE TABLE IF NOT EXISTS metabolite_enzymes (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    chemical_id   TEXT NOT NULL REFERENCES chemical_entities,
    gene_symbol   TEXT NOT NULL,
    entity_id     TEXT REFERENCES entities,
    role          TEXT NOT NULL,   -- synthesis|precursor_synthesis|degradation|efflux|influx
    source_id     TEXT NOT NULL,
    raw_role      TEXT,
    direction_known INTEGER DEFAULT 1,
    UNIQUE (chemical_id, gene_symbol, role, source_id)
);

CREATE TABLE IF NOT EXISTS chem_source_records (
    record_id           INTEGER PRIMARY KEY AUTOINCREMENT,
    source_id           TEXT NOT NULL,
    chem_interaction_id TEXT,
    chemical_id         TEXT,
    raw_ligand          TEXT,      -- verbatim, INCLUDING composite strings
    raw_receptor        TEXT,
    raw_reference       TEXT,
    raw_row             TEXT       -- JSON of the original row
);

CREATE INDEX IF NOT EXISTS ix_chemsrc_inter ON chem_source_records(chem_interaction_id);

-- Work table produced by np3 and consumed by np4. One row per
-- (source, raw receptor token, resolved entity). A receptor complex explodes
-- into one row per subunit, tied together by complex_id (decision 6).
CREATE TABLE IF NOT EXISTS np_receptor_map (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    source_id      TEXT NOT NULL,
    raw_receptor   TEXT NOT NULL,
    entity_id      TEXT REFERENCES entities,
    symbol         TEXT,
    uniprot        TEXT,
    complex_id     TEXT,
    complex_role   TEXT,
    receptor_type  TEXT,
    resolution     TEXT,          -- matched_symbol | matched_uniprot | created | unresolved
    UNIQUE (source_id, raw_receptor, entity_id)
);

-- Evidence for the chemical half. Same shape as the peptide `evidence` table and
-- the same rules: one row per (pair, paper); a review is a lead, not a receipt;
-- nothing is 'confirmed' without a verbatim sentence naming both partners.
-- Kept SEPARATE from `evidence` so the two halves can never be scored together.
CREATE TABLE IF NOT EXISTS np_evidence (
    np_evidence_id      TEXT PRIMARY KEY,
    chem_interaction_id TEXT NOT NULL REFERENCES chemical_interactions,
    pmid                TEXT NOT NULL,
    source_id           TEXT NOT NULL,
    verification_outcome TEXT NOT NULL DEFAULT 'pending',
    evidence_depth      TEXT,
    evidence_sentence   TEXT,
    assay_class         TEXT,
    effect_direction    TEXT,
    extraction_method   TEXT,
    extraction_confidence TEXT,
    notes               TEXT,
    verified_at         TEXT,
    UNIQUE (chem_interaction_id, pmid)
);
CREATE INDEX IF NOT EXISTS ix_npev_pmid ON np_evidence(pmid);
CREATE INDEX IF NOT EXISTS ix_npev_outcome ON np_evidence(verification_outcome);
