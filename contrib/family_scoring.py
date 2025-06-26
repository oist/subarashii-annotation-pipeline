#!/usr/bin/env python3
"""
family_scoring.py  –  rank trimmed alignments for concatenation markers
and SpeciesRax gene-tree inference.

Usage:
    python family_scoring.py  [inflation] [family_set]
Defaults:
    inflation   = 1.8
    family_set  = Normal
Output:
    results/families/<inflation>/family_scores.tsv  (tab-separated)
Columns:
    family, n_species, n_total, median_len, single_copy,
    dist_to_2, oversize_penalty,
    concat_score, speciesrax_score
"""
import sys, subprocess, pathlib, statistics, math
from collections import Counter, defaultdict
from Bio import SeqIO

# -------- CLI ----------------------------------------------------
infl   = sys.argv[1] if len(sys.argv) > 1 else "1.8"
subset = sys.argv[2] if len(sys.argv) > 2 else "Normal"

# -------- read result_dir from pipeline.conf ---------------------
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
read_conf  = SCRIPT_DIR / "read_conf.sh"
conf_file  = SCRIPT_DIR / "pipeline.conf"

def conf_get(section, key):
    cmd = ["bash", "-c",
           f"source {read_conf} && conf_get {conf_file} {section} {key}"]
    return subprocess.check_output(cmd, text=True).strip()

RESULT_DIR = pathlib.Path(conf_get("dirs", "result_dir"))
ALIGN_DIR  = RESULT_DIR / f"families/{infl}/{subset}"
OUT_TSV    = RESULT_DIR / f"families/{infl}/family_scores.tsv"

if not ALIGN_DIR.exists():
    sys.exit(f"[error] folder not found: {ALIGN_DIR}")

trim_paths = sorted(ALIGN_DIR.glob("*.trim"))
if not trim_paths:
    sys.exit(f"[error] no .trim files in {ALIGN_DIR}")

# -------- first pass: collect species universe & per-family data -
total_species_set = set()
fam_data = {}

def ungapped(rec):
    return len(str(rec.seq).replace("-", ""))

for trim in trim_paths:
    fam = trim.stem.split(".")[0]          # gf123 from gf123.trim
    recs = list(SeqIO.parse(trim, "fasta"))
    if not recs:
        continue

    counts = Counter(r.id.split("_g")[0] for r in recs)
    total_species_set.update(counts.keys())

    lengths = [ungapped(r) for r in recs]

    fam_data[fam] = dict(
        n_species=len(counts),
        n_total=len(recs),
        median_len=int(statistics.median(lengths)),
        single_copy=1 if all(v == 1 for v in counts.values()) else 0,
        dist_to_2=math.sqrt(sum((v - 2) ** 2 for v in counts.values())),
        counts=counts                          # keep for oversize calc
    )

TOTAL_SPECIES = len(total_species_set)
if TOTAL_SPECIES == 0:
    sys.exit("[error] could not detect any species codes")

# -------- scoring ------------------------------------------------
rows = []
for fam, d in fam_data.items():
    oversize_penalty = max(0, d["n_total"] - 3 * d["n_species"]) / d["n_species"]
    universality_fraction = d["n_species"] / TOTAL_SPECIES

    concat_score = d["single_copy"] * (10 * d["n_species"] + 0.01 * d["median_len"])

    speciesrax_score = (
        universality_fraction * (0.01 * d["median_len"] + d["n_species"])
    ) / (1 + d["dist_to_2"] + oversize_penalty)

    rows.append((
        fam,
        d["n_species"], d["n_total"], d["median_len"],
        d["single_copy"],
        round(d["dist_to_2"], 2),
        round(oversize_penalty, 2),
        round(concat_score, 2),
        round(speciesrax_score, 4)
    ))

# -------- write table -------------------------------------------
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
with OUT_TSV.open("w") as f:
    f.write("\t".join([
        "family", "n_species", "n_total", "median_len",
        "single_copy", "dist_to_2", "oversize_penalty",
        "concat_score", "speciesrax_score"]) + "\n")
    for r in sorted(rows, key=lambda x: x[-1], reverse=True):
        f.write("\t".join(map(str, r)) + "\n")

print(f"[family_scoring] wrote {OUT_TSV} ({len(rows)} families)")

