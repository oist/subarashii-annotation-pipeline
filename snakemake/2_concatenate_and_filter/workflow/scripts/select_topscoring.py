#!/usr/bin/env python3
"""
Prepare two alignment sets:

1. Concatenation candidates
   - universality >= 0.80
   - keep .trim
   - concatenate and write concatenated.faa + concatenated.phy

2. PhyloBayes candidates
   - top 300 by speciesrax_score
   - copy .trim, convert each to sequential PHYLIP

Usage
  python 16_prepare_alignments.py [inflation] [subset]

Defaults
  inflation = 1.8
  subset    = Normal
"""

import sys, subprocess, statistics, pathlib, shutil
from collections import defaultdict
import pandas as pd
from Bio import SeqIO

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def write_phylip_relaxed(records, outfile):
    """
    records : iterable of SeqRecord
    outfile : Path | str
    Writes PHYLIP-relaxed (seq type inferred).
    """
    aln = MultipleSeqAlignment(list(records))
    AlignIO.write(aln, outfile, "phylip-relaxed")



infl   = sys.argv[1] if len(sys.argv) > 1 else "1.8"
subset = sys.argv[2] if len(sys.argv) > 2 else "Normal"

# --- locate result_dir via read_conf.sh --------------------------
SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
read_conf  = SCRIPT_DIR / "read_conf.sh"
conf_file  = SCRIPT_DIR / "pipeline.conf"

def conf_get(sec, key):
    cmd=["bash","-c",f"source {read_conf} && conf_get {conf_file} {sec} {key}"]
    return subprocess.check_output(cmd, text=True).strip()

RES = pathlib.Path(conf_get("dirs","result_dir"))
ALIGN_DIR = RES / f"families/{infl}/{subset}"
SCORE_TSV = RES / f"families/{infl}/family_scores.tsv"

CONCAT_DIR = RES / f"families/{infl}/Concatenation"
PB_DIR     = RES / f"families/{infl}/PhyloBayes"
CONCAT_DIR.mkdir(parents=True, exist_ok=True)
PB_DIR.mkdir(parents=True, exist_ok=True)

# --- load scores --------------------------------------------------
df = pd.read_csv(SCORE_TSV, sep="\t")
df["universality"] = df["n_species"] / df["n_species"].max()

# ---------- 1. Concatenation set ---------------------------------
concat_fams = df[df["universality"] >= 0.80]["family"].tolist()
print(f"[concat] {len(concat_fams)} families ≥80 % universality")

# first pass: discover universe of species
all_species = set()
for fam in concat_fams:
    trim = ALIGN_DIR / f"{fam}.trim"
    if not trim.exists():
        continue
    for r in SeqIO.parse(trim, "fasta"):
        all_species.add(r.id.split("_g")[0])

all_species = sorted(all_species)
concatenated_seqs = {sp: "" for sp in all_species}

for fam in concat_fams:
    trim = ALIGN_DIR / f"{fam}.trim"
    if not trim.exists():
        continue

    # copy .trim once
    dst = CONCAT_DIR / trim.name
    if not dst.exists():
        shutil.copy(trim, dst)

    recs = list(SeqIO.parse(trim, "fasta"))
    fam_len = len(recs[0].seq)
    seq_by_species = {r.id.split("_g")[0]: str(r.seq) for r in recs}

    for sp in all_species:
        concatenated_seqs[sp] += seq_by_species.get(sp, "-" * fam_len)

# ------------- write concatenated fasta --------------------------
concat_faa = CONCAT_DIR / "concatenated.faa"
with concat_faa.open("w") as f:
    for sp in all_species:
        f.write(f">{sp}\n{concatenated_seqs[sp]}\n")
print(f"[concat] wrote {concat_faa}")

# ------------- write concatenated PHYLIP-relaxed -----------------
concat_phy = CONCAT_DIR / "concatenated.phy"
records = [SeqRecord(Seq(seq), id=sp, description="") 
           for sp, seq in concatenated_seqs.items()]
write_phylip_relaxed(records, concat_phy)
print(f"[concat] wrote {concat_phy} (PHYLIP-relaxed)")




# ---------- 2. PhyloBayes set (top 300) --------------------------
top_pb = df.sort_values("speciesrax_score", ascending=False).head(300)

for fam in top_pb["family"]:
    trim = ALIGN_DIR / f"{fam}.trim"
    if not trim.exists(): continue
    dst_dir = PB_DIR / fam
    dst_dir.mkdir(exist_ok=True)
    dst_trim = dst_dir / trim.name
    if not dst_trim.exists(): shutil.copy(trim, dst_trim)

    # convert to sequential phylip
    phy_path = dst_dir / f"{fam}.phy"
    recs = list(SeqIO.parse(trim, "fasta"))
    write_phylip_relaxed(recs, phy_path)

print(f"[pb] prepared {len(top_pb)} families under {PB_DIR}")

