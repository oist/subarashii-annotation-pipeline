#!/usr/bin/env python3
"""
Create a DuckDB view layer + Parquet files from eggnog-mapper pipeline output.

Outputs:
  <output_dir>/genome.parquet      — one row per genome (taxonomy + GTDB quality metrics)
  <output_dir>/protein.parquet     — one row per protein (eggnog OGs, KEGG, aa_seq,
                                     genomic coordinates)
  <output_dir>/protein_dna.parquet — one row per protein (CDS nucleotide sequence),
                                     kept separate because DNA seqs are ~3× larger
  <duckdb>                         — DuckDB with genome table, protein view,
                                     and protein_dna view

The CDS nucleotide sequences are extracted from the selected nucleotide genome FASTA
files using the coordinates stored in the Prodigal protein headers
(sequenceid_new2original.csv). Reverse-complement is applied for minus-strand genes.

Usage:
  scripts/create_db.py \
      --genome2abbrev  results/{dataset}/genome2abbrev.csv \
      --eggnog         results/{dataset}/eggnog/annotations.emapper.annotations \
      --fasta-dir      results/{dataset}/genomes_protein \
      --na-dir         resources/{dataset}/selected_genomes_na \
      --id-map         results/{dataset}/sequenceid_new2original.csv \
      --bac-metadata   resources/bac120_metadata_r226.tsv.gz \
      --arc-metadata   resources/ar53_metadata_r226.tsv.gz \
      --output-dir     results/{dataset}/db \
      --duckdb         results/{dataset}/db/dataset.duckdb

Authors:
    - Claude Code Sonnet

Version:
    - v0.1 (2026-09-02): initial concept
"""

import argparse
import csv
import gzip
import os
import re
import sys
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO


# ── taxonomy ───────────────────────────────────────────────────────────────────

TAXONOMY_PREFIXES = [
    ("d__", "domain"),
    ("p__", "phylum"),
    ("c__", "class"),
    ("o__", "order"),
    ("f__", "family"),
    ("g__", "genus"),
    ("s__", "species"),
]


def parse_taxonomy(tax_str: str) -> dict:
    """Parse 'd__Bacteria;p__...;s__...' into a dict keyed by rank name."""
    result = {col: "" for _, col in TAXONOMY_PREFIXES}
    for part in tax_str.split(";"):
        part = part.strip()
        for prefix, col in TAXONOMY_PREFIXES:
            if part.startswith(prefix):
                result[col] = part[len(prefix):]
    return result


# ── GTDB metadata ──────────────────────────────────────────────────────────────

# Columns we want from the GTDB metadata TSV, with per-version name fallbacks.
# First match wins.
METADATA_COLUMN_MAP = {
    "completeness":      ["checkm2_completeness",  "checkm_completeness"],
    "contamination":     ["checkm2_contamination",  "checkm_contamination"],
    "coding_density":    ["coding_density"],
    "contig_count":      ["contig_count"],
    "gc_percentage":     ["gc_percentage"],
    "genome_size":       ["genome_size"],
    "gtdb_representative": ["gtdb_representative"],
    "n50_contigs":       ["n50_contigs"],
    "ncbi_organism_name": ["ncbi_organism_name"],
    "protein_count":     ["protein_count"],
    "ssu_count":         ["ssu_count"],
    "trna_count":        ["trna_count", "trna_aa_count"],
}

METADATA_FLOAT_COLS  = {"completeness", "contamination", "coding_density", "gc_percentage"}
METADATA_INT_COLS    = {"contig_count", "genome_size", "n50_contigs", "protein_count",
                         "ssu_count", "trna_count"}
METADATA_BOOL_COLS   = {"gtdb_representative"}


def _open_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_gtdb_metadata(bac_path: str, arc_path: str, acc_set: set) -> dict:
    """
    Load GTDB quality metrics for the genomes in acc_set from the metadata TSVs.

    Returns {accession: {field: value, ...}}.  Missing values become None.
    """
    result = {}
    for path in (bac_path, arc_path):
        if not path or not os.path.exists(path):
            print(f"[db] warning: metadata file not found: {path}", file=sys.stderr)
            continue
        with _open_maybe_gz(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            headers = reader.fieldnames or []
            # resolve which source column maps to each output field
            col_map = {}
            for out_col, candidates in METADATA_COLUMN_MAP.items():
                for c in candidates:
                    if c in headers:
                        col_map[out_col] = c
                        break
            for row in reader:
                acc = row.get("accession", "").strip()
                if acc not in acc_set:
                    continue
                rec = {}
                for out_col, src_col in col_map.items():
                    raw = row.get(src_col, "").strip()
                    if raw in ("", "N/A", "none", "None", "na", "NA"):
                        rec[out_col] = None
                    elif out_col in METADATA_FLOAT_COLS:
                        try:
                            rec[out_col] = float(raw)
                        except ValueError:
                            rec[out_col] = None
                    elif out_col in METADATA_INT_COLS:
                        try:
                            rec[out_col] = int(float(raw))
                        except ValueError:
                            rec[out_col] = None
                    elif out_col in METADATA_BOOL_COLS:
                        rec[out_col] = raw.lower() in ("true", "t", "yes", "1")
                    else:
                        rec[out_col] = raw or None
                result[acc] = rec
    return result


# ── genome2abbrev ──────────────────────────────────────────────────────────────

def load_genome_abbrev(path: str) -> tuple[list[dict], dict, dict]:
    """
    Load genome2abbrev.csv (header: accession,short,taxa).

    Returns (genome_rows_base, acc_to_short, short_to_acc).
    """
    genome_rows = []
    acc_to_short = {}
    short_to_acc = {}

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            acc   = row["accession"]
            short = row["short"]
            taxa  = parse_taxonomy(row.get("taxa", ""))
            acc_to_short[acc]   = short
            short_to_acc[short] = acc
            genome_rows.append({"acc": acc, "short": short, **taxa})

    return genome_rows, acc_to_short, short_to_acc


GENOME_SCHEMA = pa.schema([
    pa.field("acc",              pa.string()),
    pa.field("short",            pa.string()),
    pa.field("domain",           pa.string()),
    pa.field("phylum",           pa.string()),
    pa.field("class",            pa.string()),
    pa.field("order",            pa.string()),
    pa.field("family",           pa.string()),
    pa.field("genus",            pa.string()),
    pa.field("species",          pa.string()),
    pa.field("completeness",     pa.float64()),
    pa.field("contamination",    pa.float64()),
    pa.field("coding_density",   pa.float64()),
    pa.field("contig_count",     pa.int64()),
    pa.field("gc_percentage",    pa.float64()),
    pa.field("genome_size",      pa.int64()),
    pa.field("gtdb_representative", pa.bool_()),
    pa.field("n50_contigs",      pa.int64()),
    pa.field("ncbi_organism_name", pa.string()),
    pa.field("protein_count",    pa.int64()),
    pa.field("ssu_count",        pa.int64()),
    pa.field("trna_count",       pa.int64()),
])

EMPTY_META = {c: None for c in [
    "completeness", "contamination", "coding_density", "contig_count",
    "gc_percentage", "genome_size", "gtdb_representative", "n50_contigs",
    "ncbi_organism_name", "protein_count", "ssu_count", "trna_count",
]}


def write_genome_parquet(genome_rows: list[dict], metadata: dict, out_path: str) -> None:
    rows = []
    for g in genome_rows:
        meta = metadata.get(g["acc"], EMPTY_META)
        rows.append({**g, **meta})
    table = pa.Table.from_pylist(rows, schema=GENOME_SCHEMA)
    pq.write_table(table, out_path)
    n_with_meta = sum(1 for g in genome_rows if g["acc"] in metadata)
    print(f"[db] genome.parquet: {len(rows)} genomes "
          f"({n_with_meta} with GTDB quality metrics)")


# ── protein ID mapping (sequenceid_new2original.csv) ──────────────────────────

_PRODIGAL_RE = re.compile(
    r"(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>-?\d+)\s+#\s+(?P<attrs>[^#]*)"
)
_PARTIAL_RE = re.compile(r"partial=(\d{2})")


def parse_protein_header(description: str) -> dict:
    """
    Parse a GTDB/Prodigal FASTA header description into location fields.

    Header format (after the protein ID):
      {protein_id} # {start} # {end} # {strand} # ID=...;partial=NN;...

    Returns {"contig": str, "start": int, "end": int, "strand": int, "partial": str}.
    All fields are None on parse failure.
    """
    # Split on first ' # ' to separate protein_id from coordinates
    parts = description.split(" # ", 1)
    protein_part = parts[0].strip()

    # contig: remove any leading "{accession}_prot_" prefix, then strip last _N
    if "_prot_" in protein_part:
        protein_part = protein_part.split("_prot_", 1)[1]
    m_contig = re.match(r"^(.+)_\d+$", protein_part)
    contig = m_contig.group(1) if m_contig else None

    if len(parts) < 2:
        return {"contig": contig, "start": None, "end": None,
                "strand": None, "partial": None}

    m = _PRODIGAL_RE.search(parts[1])
    if not m:
        return {"contig": contig, "start": None, "end": None,
                "strand": None, "partial": None}

    mp = _PARTIAL_RE.search(m.group("attrs"))
    return {
        "contig":  contig,
        "start":   int(m.group("start")),
        "end":     int(m.group("end")),
        "strand":  int(m.group("strand")),
        "partial": mp.group(1) if mp else None,
    }


def load_id_map(path: str) -> dict:
    """
    Load sequenceid_new2original.csv (actually TSV: new_id \\t old_description).

    Returns {genome_short: {protein_id: location_dict}}.

    Grouping by genome short-code avoids an O(N_genomes × N_proteins) scan when
    looking up proteins during per-genome iteration.
    """
    _short_re = re.compile(r"^(.+)_g\d+$")
    id_map: dict[str, dict] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            new_id, _, description = line.partition("\t")
            new_id = new_id.strip()
            m = _short_re.match(new_id)
            short = m.group(1) if m else new_id
            id_map.setdefault(short, {})[new_id] = parse_protein_header(description)
    return id_map


# ── eggnog annotations ─────────────────────────────────────────────────────────

def load_eggnog_annotations(path: str) -> dict:
    """
    Parse emapper.annotations into {protein_id: annotation_dict}.

    Only the fields we write to the parquet are kept.
    """
    annotations = {}
    col = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                headers = line.lstrip("#").rstrip("\n").split("\t")
                col = {h.strip(): i for i, h in enumerate(headers)}
                continue
            if not col:
                continue

            parts = line.rstrip("\n").split("\t")
            if not parts or parts[0].startswith("#"):
                continue

            def get(name, default=None):
                i = col.get(name)
                if i is None or i >= len(parts):
                    return default
                v = parts[i].strip()
                return None if v in ("", "-") else v

            protein_id = get("query") or get("#query") or parts[0].strip()

            ogs_str = get("eggNOG_OGs") or ""
            ogs_list = [og.strip() for og in ogs_str.split(",") if og.strip()]

            annotations[protein_id] = {
                "eggnog7_ogs":    ogs_list,
                "cog_category":   get("COG_category"),
                "description":    get("Description"),
                "preferred_name": get("Preferred_name"),
                "kegg_ko":        get("KEGG_ko"),
                "go":             get("GOs"),
                "pfam":           get("PFAMs"),
            }

    return annotations


# ── protein parquet (streaming per genome) ────────────────────────────────────

PROTEIN_SCHEMA = pa.schema([
    pa.field("protein_id",     pa.string()),
    pa.field("acc",            pa.string()),
    pa.field("eggnog7_ogs",    pa.list_(pa.string())),
    pa.field("cog_category",   pa.string()),
    pa.field("description",    pa.string()),
    pa.field("preferred_name", pa.string()),
    pa.field("kegg_ko",        pa.string()),
    pa.field("go",             pa.string()),
    pa.field("pfam",           pa.string()),
    pa.field("aa_seq",         pa.string()),
    pa.field("contig",         pa.string()),
    pa.field("start",          pa.int32()),
    pa.field("end",            pa.int32()),
    pa.field("strand",         pa.int8()),
    pa.field("partial",        pa.string()),
])

_EMPTY_ANNOT = {
    "eggnog7_ogs": [], "cog_category": None, "description": None,
    "preferred_name": None, "kegg_ko": None, "go": None, "pfam": None,
}
_EMPTY_LOC = {
    "contig": None, "start": None, "end": None, "strand": None, "partial": None,
}


def write_protein_parquet(
    short_to_acc: dict,
    fasta_dir: str,
    annotations: dict,
    id_map: dict,
    out_path: str,
    row_group_size: int,
    compression: str,
) -> int:
    """
    Stream protein records genome by genome into a parquet file.

    id_map must be the grouped structure returned by load_id_map:
        {short: {protein_id: loc_dict}}

    row_group_size and compression are passed through to ParquetWriter so
    the SLURM memory budget drives how much data is buffered before each flush.
    """
    total = 0
    writer = pq.ParquetWriter(out_path, PROTEIN_SCHEMA, compression=compression)

    faa_files = sorted(Path(fasta_dir).glob("*.faa"))
    for faa_path in faa_files:
        short       = faa_path.stem
        acc         = short_to_acc.get(short, short)
        genome_locs = id_map.get(short, {})
        rows        = []

        for rec in SeqIO.parse(str(faa_path), "fasta"):
            pid   = rec.id
            annot = annotations.get(pid, _EMPTY_ANNOT)
            loc   = genome_locs.get(pid, _EMPTY_LOC)

            rows.append({
                "protein_id":     pid,
                "acc":            acc,
                "eggnog7_ogs":    annot["eggnog7_ogs"],
                "cog_category":   annot["cog_category"],
                "description":    annot["description"],
                "preferred_name": annot["preferred_name"],
                "kegg_ko":        annot["kegg_ko"],
                "go":             annot["go"],
                "pfam":           annot["pfam"],
                "aa_seq":         str(rec.seq),
                "contig":         loc["contig"],
                "start":          loc["start"],
                "end":            loc["end"],
                "strand":         loc["strand"],
                "partial":        loc["partial"],
            })

            # Flush when the batch hits the row group budget
            if len(rows) >= row_group_size:
                writer.write_batch(
                    pa.RecordBatch.from_pylist(rows, schema=PROTEIN_SCHEMA)
                )
                total += len(rows)
                rows = []

        if rows:
            writer.write_batch(pa.RecordBatch.from_pylist(rows, schema=PROTEIN_SCHEMA))
            total += len(rows)

    writer.close()
    print(f"[db] protein.parquet: {total} proteins from {len(faa_files)} genomes")
    return total


# ── CDS nucleotide extraction ──────────────────────────────────────────────────

PROTEIN_DNA_SCHEMA = pa.schema([
    pa.field("protein_id", pa.string()),
    pa.field("acc",        pa.string()),
    pa.field("na_seq",     pa.string()),
])


def _build_contig_index(na_fasta_path: str) -> dict:
    """Return {contig_id: Bio.Seq} for all contigs in a nucleotide FASTA."""
    return {rec.id: rec.seq for rec in SeqIO.parse(na_fasta_path, "fasta")}


def write_protein_dna_parquet(
    short_to_acc: dict,
    na_dir: str,
    id_map: dict,
    out_path: str,
    row_group_size: int,
    compression: str,
) -> int:
    """
    Extract CDS nucleotide sequences from per-genome nucleotide FASTA files.

    id_map must be the grouped structure returned by load_id_map:
        {short: {protein_id: loc_dict}}

    Grouping means we only iterate each genome's own proteins — no full-dataset
    scan per genome.  Coordinates are 1-based (Prodigal); forward strand uses
    dna[start-1:end], reverse strand is reverse-complemented.
    """
    total = 0
    writer = pq.ParquetWriter(out_path, PROTEIN_DNA_SCHEMA, compression=compression)

    na_dir_path = Path(na_dir)
    na_files = sorted(
        f for f in na_dir_path.iterdir()
        if f.suffix in (".fna", ".fa", ".fasta") or f.name.endswith(".fna.gz")
    )

    for na_path in na_files:
        short = na_path.name
        for ext in (".fna.gz", ".fna", ".fa", ".fasta"):
            if short.endswith(ext):
                short = short[: -len(ext)]
                break
        short = re.sub(r"_genomic$", "", short)

        acc = short_to_acc.get(short)
        if acc is None:
            continue

        genome_locs = id_map.get(short, {})
        if not genome_locs:
            continue

        contig_index = _build_contig_index(str(na_path))
        rows = []

        for pid, loc in genome_locs.items():
            contig = loc.get("contig")
            start  = loc.get("start")
            end    = loc.get("end")
            strand = loc.get("strand")

            if contig is None or start is None or end is None:
                rows.append({"protein_id": pid, "acc": acc, "na_seq": None})
            else:
                contig_seq = contig_index.get(contig)
                if contig_seq is None:
                    rows.append({"protein_id": pid, "acc": acc, "na_seq": None})
                else:
                    cds = contig_seq[start - 1: end]   # 1-based → 0-based
                    if strand == -1:
                        cds = cds.reverse_complement()
                    rows.append({"protein_id": pid, "acc": acc, "na_seq": str(cds)})

            if len(rows) >= row_group_size:
                writer.write_batch(
                    pa.RecordBatch.from_pylist(rows, schema=PROTEIN_DNA_SCHEMA)
                )
                total += len(rows)
                rows = []

        if rows:
            writer.write_batch(pa.RecordBatch.from_pylist(rows, schema=PROTEIN_DNA_SCHEMA))
            total += len(rows)

    writer.close()
    print(f"[db] protein_dna.parquet: {total} CDS sequences from {len(na_files)} genomes")
    return total


# ── DuckDB ─────────────────────────────────────────────────────────────────────

def create_duckdb(
    db_path: str,
    genome_parquet: str,
    protein_parquet: str,
    protein_dna_parquet: str,
) -> None:
    genome_abs      = os.path.abspath(genome_parquet)
    protein_abs     = os.path.abspath(protein_parquet)
    protein_dna_abs = os.path.abspath(protein_dna_parquet)

    if os.path.exists(db_path):
        os.remove(db_path)

    con = duckdb.connect(db_path)
    con.execute(f"CREATE TABLE genome AS SELECT * FROM read_parquet('{genome_abs}')")
    con.execute(f"CREATE VIEW protein     AS SELECT * FROM read_parquet('{protein_abs}')")
    con.execute(f"CREATE VIEW protein_dna AS SELECT * FROM read_parquet('{protein_dna_abs}')")
    con.close()

    def mb(p): return os.path.getsize(p) // 1024 // 1024

    print(f"[db] DuckDB: {db_path}")
    print(f"[db]   genome      — base table  ({os.path.getsize(genome_parquet)//1024} KB)")
    print(f"[db]   protein     — view  ({mb(protein_parquet)} MB)")
    print(f"[db]   protein_dna — view  ({mb(protein_dna_parquet)} MB)")
    print(f"[db]   Example queries:")
    print(f"[db]     SELECT * EXCLUDE (aa_seq) FROM protein LIMIT 5;")
    print(f"[db]     -- join AA and CDS sequences:")
    print(f"[db]     SELECT p.protein_id, p.aa_seq, d.na_seq")
    print(f"[db]       FROM protein p JOIN protein_dna d USING (protein_id) LIMIT 5;")


# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--genome2abbrev",  required=True,
                        help="genome2abbrev.csv (accession,short,taxa)")
    parser.add_argument("--eggnog",         required=True,
                        help="emapper.annotations file")
    parser.add_argument("--fasta-dir",      required=True,
                        help="Directory of per-genome renamed AA FASTA files ({short}.faa)")
    parser.add_argument("--na-dir",         required=True,
                        help="Directory of per-genome nucleotide FASTA files "
                             "(resources/{dataset}/selected_genomes_na)")
    parser.add_argument("--id-map",         required=True,
                        help="sequenceid_new2original.csv (new_id TAB original_header)")
    parser.add_argument("--bac-metadata",   required=True,
                        help="GTDB bacterial metadata TSV (or .tsv.gz)")
    parser.add_argument("--arc-metadata",   required=True,
                        help="GTDB archaeal metadata TSV (or .tsv.gz)")
    parser.add_argument("--output-dir",     required=True,
                        help="Directory for parquet files")
    parser.add_argument("--duckdb",         required=True,
                        help="Path for the output .duckdb file")
    parser.add_argument("--max-ram-mb",     type=int, default=8000,
                        help="RAM budget in MB (pass SLURM mem_mb_per_cpu here). "
                             "Controls ParquetWriter row-group sizes so each flush "
                             "stays within ~20%% of this budget. Default: 8000.")
    parser.add_argument("--compression",    default="zstd",
                        choices=["zstd", "snappy", "gzip", "brotli", "none"],
                        help="Parquet compression codec. zstd (default) gives the "
                             "best size/speed trade-off; snappy is faster but larger.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Row-group size: target 20%% of the RAM budget per flush.
    # Estimated uncompressed bytes per row:
    #   protein.parquet     — aa_seq ~300 B + eggnog/other ~200 B  ≈ 500 B
    #   protein_dna.parquet — na_seq ~900 B + ids ~50 B            ≈ 950 B
    budget_bytes   = args.max_ram_mb * 1024 * 1024 * 0.20
    rg_protein     = max(1_000, int(budget_bytes / 500))
    rg_protein_dna = max(1_000, int(budget_bytes / 950))
    print(f"[db] RAM budget: {args.max_ram_mb} MB  "
          f"→ row-group sizes: protein={rg_protein:,}  protein_dna={rg_protein_dna:,}  "
          f"compression={args.compression}")

    # 1. Genome abbreviation map + base taxonomy
    print("[db] Loading genome abbreviation map …")
    genome_rows, acc_to_short, short_to_acc = load_genome_abbrev(args.genome2abbrev)
    acc_set = set(g["acc"] for g in genome_rows)

    # 2. GTDB quality metrics
    print("[db] Loading GTDB metadata …")
    metadata = load_gtdb_metadata(args.bac_metadata, args.arc_metadata, acc_set)

    genome_parquet = os.path.join(args.output_dir, "genome.parquet")
    write_genome_parquet(genome_rows, metadata, genome_parquet)

    # 3. eggnog annotations
    print("[db] Loading eggnog-mapper annotations …")
    annotations = load_eggnog_annotations(args.eggnog)
    print(f"[db]   {len(annotations)} annotated proteins")

    # 4. Protein ID → location map (grouped by genome short code)
    print("[db] Loading protein ID map …")
    id_map = load_id_map(args.id_map)
    n_proteins = sum(len(v) for v in id_map.values())
    print(f"[db]   {n_proteins} protein ID mappings across {len(id_map)} genomes")

    # 5. Protein parquet with AA sequences (streamed per genome)
    print("[db] Writing protein.parquet (includes aa_seq — may be large) …")
    protein_parquet = os.path.join(args.output_dir, "protein.parquet")
    write_protein_parquet(
        short_to_acc, args.fasta_dir, annotations, id_map,
        protein_parquet, rg_protein, args.compression,
    )

    # 6. Protein DNA parquet — CDS extracted from nucleotide FASTA by coordinates
    print("[db] Writing protein_dna.parquet (CDS sequences — will be larger) …")
    protein_dna_parquet = os.path.join(args.output_dir, "protein_dna.parquet")
    write_protein_dna_parquet(
        short_to_acc, args.na_dir, id_map,
        protein_dna_parquet, rg_protein_dna, args.compression,
    )

    # 7. DuckDB view layer
    print("[db] Creating DuckDB …")
    create_duckdb(args.duckdb, genome_parquet, protein_parquet, protein_dna_parquet)

    print("[db] Done.")


if __name__ == "__main__":
    main()
