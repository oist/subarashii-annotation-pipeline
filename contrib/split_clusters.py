#!/usr/bin/env python3
"""
Split one MCL cluster file into separate family FASTA files.

Reads settings from pipeline.conf / pipeline.local.conf:
  [dirs]   result_dir
  [families] inflation  (cluster file comes from mcl/clusters_I<inflation>.txt)
              huge_cutoff

Creates subfolders in <result_dir>/families/<inflation>/:
   Singletons / Duplets / Triplets / Normal / Huge

Also writes two logs in the same directory:
   family_type_counts.tsv
   family_details.tsv
"""

import configparser
import pathlib
import statistics
from collections import Counter, defaultdict
from Bio import SeqIO

# ------------------------- read config -------------------------
cfg = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))
cfg.read(["pipeline.conf", "pipeline.local.conf"])

res_dir     = pathlib.Path(cfg["dirs"]["result_dir"])
infl        = cfg["families"].get("inflation", "1.8")
huge_cutoff = int(cfg["families"].get("huge_cutoff", "300"))

cluster_file = res_dir / "mcl" / f"clusters_I{infl}.txt"
out_root     = res_dir / "families" / infl
out_root.mkdir(parents=True, exist_ok=True)

# ------------------------- load sequences once -----------------
fasta_path = res_dir / "all_proteins.faa"
seq_dict   = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}
print(f"[split_clusters] loaded {len(seq_dict)} sequences")

# ------------------------- helpers -----------------------------
def family_folder(size):
    if size == 1:
        return "Singletons"
    if size == 2:
        return "Duplets"
    if size == 3:
        return "Triplets"
    if size >= huge_cutoff:
        return "Huge"
    return "Normal"

def lengths(records):
    return [len(r.seq) for r in records]

def has_paralog(records):
    species_seen = Counter(r.id.split("_g")[0] for r in records)
    return any(cnt > 1 for cnt in species_seen.values())

# ------------------------- iterate clusters --------------------
type_counter = Counter()
detail_rows  = []

with cluster_file.open() as fh:
    for idx, line in enumerate(fh):
        ids = line.strip().split()
        size = len(ids)
        folder = family_folder(size)
        type_counter[folder] += 1

        gf_name = f"gf{idx}"
        out_dir = out_root / folder
        out_dir.mkdir(exist_ok=True)
        out_faa = out_dir / f"{gf_name}.faa"

        recs = [seq_dict[i] for i in ids if i in seq_dict]
        SeqIO.write(recs, out_faa.open("w"), "fasta")

        lens = lengths(recs)
        detail_rows.append((
            gf_name, size,
            max(lens), min(lens),
            int(statistics.median(lens)),
            "Yes" if has_paralog(recs) else "No"
        ))

# ------------------------- write logs --------------------------
with (out_root / "family_type_counts.tsv").open("w") as f:
    for t in ("Singletons", "Duplets", "Triplets", "Normal", "Huge"):
        f.write(f"{t}\t{type_counter[t]}\n")

with (out_root / "family_details.tsv").open("w") as f:
    f.write("family\tsize\tlongest\tshortest\tmedian\tparalog_in_species\n")
    for row in detail_rows:
        f.write("\t".join(map(str, row)) + "\n")

print("[split_clusters] done. Outputs in", out_root)

