# EggNOG-mapper v3 Output Format

Notes on how EggNOG-mapper v3 (beta) output differs from v2, what each field
means, and how it is used by this pipeline.

## Column changes relative to v2

EggNOG-mapper v3 introduced a new column layout.  The two columns that changed
most significantly for this pipeline are:

| Column (0-based) | Header | v2 content | v3 content |
|---|---|---|---|
| 4 | `eggNOG_OGs` | COG IDs at root level, e.g. `COG0457@1\|root` | Domain-family OGs, e.g. `TatC@131567\|A-1` |
| 7 | `COG_category` | Single-letter category code, e.g. `J` | Full COG family ID, e.g. `COG0805`, or a letter if no COG assigned |

The pipeline scripts (`create_cog_clusters.py`, `generate_presence_absence_table.py`)
detect the format automatically by reading the header and falling back from the
v2 strategy to the v3 strategy when no COG IDs are found in column 4.

---

## The `eggNOG_OGs` column (column 4)

### Format

```
DomainFamily@TaxID|ClusterLabel[*][,DomainFamily@TaxID|ClusterLabel ...]
```

A gene may have multiple comma-separated OG entries, one per taxonomic scope.

### Example (from a real v3 annotation file)

```
TatC@131567|A-1
B12-binding@131567|BiL-30
CsgG@131567|A-1,CsgG@1843491|Ks-13
MotA_ExbB@909932|BNU-29,MotA_ExbB@131567|RT-20
Lactamase_B@131567|HQ-13,Flavodoxin_1@131567|B-2
```

### What each part means

#### `DomainFamily` — prefix before `@`

The name of the Pfam / HMM domain family used to define this orthologous group
(e.g. `TatC`, `Lactamase_B`, `Ribosomal_S18`).  In v2 this was always a COG
identifier; in v3 it is a domain family name.

Genes with no recognisable Pfam domain use the prefix `UNK.` followed by a
short hash derived from the cluster's representative sequence, e.g.
`UNK.8IWM@28221|A-1`.  The hash is stable within a database version but has no
biological meaning beyond identifying the cluster.

#### `@TaxID` — NCBI taxonomy node

The taxonomic scope at which this orthologous group was defined.  A lower
(broader) node means the OG is more widely conserved.

Verified via NCBI Taxonomy eutils (`efetch.fcgi?db=taxonomy&id=...`):

| TaxID | Name | Rank | Notes |
|---|---|---|---|
| 1 | root | — | includes viruses |
| **131567** | **cellular organisms** | no formal rank | all Bacteria + Archaea + Eukaryota; excludes viruses |
| 2 | Bacteria | domain | |
| 2157 | Archaea | domain | |
| 2759 | Eukaryota | domain | |
| 28221 | Deltaproteobacteria | class | example of a class-level OG |
| 909932 | Negativicutes | class | example of a class-level OG |
| 1843491 | Selenomonadaceae | family | example of a family-level OG |

**Rule of thumb for phylogenomics:**
- `@131567` → universally conserved across cellular life — good universal marker gene candidates
- `@2` / `@2157` → domain-specific (bacteria-only or archaea-only)
- Narrower TaxIDs (high numbers, species/family level) → lineage-specific expansions

A gene with an OG at a very narrow taxonomic scope is less useful as a
phylogenetic marker but may be interesting for functional or comparative studies.

#### `|ClusterLabel` — internal cluster ID

The label after `|` (e.g. `A-1`, `BiL-30`, `Ks-13`, `RT-20`) is an **internal
cluster identifier** stored in the eggNOG database for a specific OG within that
domain family at that taxonomic node.

**Are the labels stable?**  
Experimentally confirmed: running the same input twice against the same
`eggnog.db` produces identical labels.  The labels are looked up from the
pre-built database, not re-computed per run, so they are:

- **Deterministic** — identical output for the same input + same database version.
- **Comparable across runs** — the same label in two annotation files (same db
  version) refers to the same reference OG and can be used as a shared key.
- **Not guaranteed stable across database versions** — if eggNOG releases a new
  database (e.g. eggNOG v6 → v7), OGs may be split, merged, or renumbered.
  Always record which `eggnog.db` version was used (visible in the `##` header
  lines of the annotation file).

**What the labels do not encode:**
- `A-1` does not mean "Archaea group 1" or any taxonomic name.
- The letter/number scheme is opaque — it carries no biological semantics beyond
  uniquely identifying a cluster within `(DomainFamily, TaxID)`.
- A trailing `*` marks the OG as the best-scoring hit for that gene.

---

## The `COG_category` column (column 7)

In v3 this column can contain:

| Value | Meaning |
|---|---|
| `COGxxxx` (e.g. `COG0805`) | COG family identifier — used by the pipeline for clustering |
| Single letter (e.g. `S`, `J`, `C`) | COG functional category but no specific COG assigned |
| `-` | No COG annotation |

The pipeline uses this column to assign genes to COG families when the v2-style
`eggNOG_OGs` column does not contain a COG identifier.

---

## Using `eggNOG_OGs` to filter by conservation breadth

If you want to select only broadly conserved genes (e.g. for universal marker
gene analysis), filter rows where any OG entry has `@131567`:

```python
import csv

with open("annotations.emapper.annotations") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if row[0].startswith("#"):
            continue
        ogs = row[4]  # eggNOG_OGs column
        is_universal = any("@131567" in og for og in ogs.split(","))
```

To resolve any TaxID to its name programmatically:

```python
import urllib.request, xml.etree.ElementTree as ET

def taxid_name(taxid: int) -> str:
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=taxonomy&id={taxid}&retmode=xml")
    with urllib.request.urlopen(url) as r:
        tree = ET.parse(r)
    return tree.find(".//ScientificName").text
```

---

## Summary

| Field | Interpretable? | Useful for |
|---|---|---|
| `DomainFamily` (before `@`) | Yes | Functional annotation, domain-level grouping |
| `@TaxID` | Yes — look up in NCBI Taxonomy | Filtering by conservation breadth |
| `\|ClusterLabel` | No — internal ID only | Identity check (same label = same OG) |
| `COG_category` col 7 | Yes (if starts with COG) | COG-based clustering (preferred in this pipeline) |
