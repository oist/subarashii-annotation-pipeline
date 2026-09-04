#!/usr/bin/env python3
"""
subarashii-annotation-pipeline governor.

Run at any time -- it checks where each dataset is and continues from there.
State is persisted between invocations, so re-running after a failure resumes
from the failed step without redoing completed work.

Usage examples:
  # Start a new dataset
  ./run_pipeline.py run --dataset my_run --ids genomes.txt

  # Resume all known datasets (e.g. after a SLURM failure)
  ./run_pipeline.py run

  # Check status without running anything
  ./run_pipeline.py status

  # Use a pre-existing species tree (skips concatenated inference step, the leaf names have to be valid GTDB genome accession IDs)
  ./run_pipeline.py run --dataset my_run --species-tree /path/to/tree.nwk

  # Override config parameters
  ./run_pipeline.py run --dataset my_run --executor local --jobs 8 --set eggnog.threads=16

Authors:
    - Claude Code Sonnet
    - Lenard L. Szantho <lenard@drenal.eu>

Version:
    - v0.5 (2026-09-04): --clustering flag makes MCL optional (default: eggnog-only)
    - v0.4 (2026-09-02): try to deduce PATH to phylobayes, iqtree3 and alerax binaries and MPI modules
    - v0.3 (2026-09-02): create db from eggnog results, parameters self-documenting, shared directory for large file storage, read in accession IDs from species tree, custom accession ID - abbreviation list
    - v0.2 (2026-08-28): eggnog mapper v3 support 
    - v0.1 (2026-08-27): initial funcionality: start snakemake, check state, continue where it left off
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent
SNAKEMAKE_DIR = ROOT / "snakemake"
STAGE1_DIR = SNAKEMAKE_DIR / "1_cluster"
STAGE2_DIR = SNAKEMAKE_DIR / "2_concatenate_and_filter"
STAGE3_DIR = SNAKEMAKE_DIR / "3_inference"

# Steps in execution order. Each is a key in the state file.
STEPS = [
    "stage1",
    "stage1_db",
    "transition_1to2",
    "stage2",
    "transition_2to3",
    "stage3_per_family",
    "stage3_concat",
    "stage3_alerax",
]

STEP_LABELS = {
    "stage1":            "Clustering (Stage 1)",
    "stage1_db":         "DuckDB/Parquet database of Eggnog Mapper results and other metadata",
    "transition_1to2":   "Transition 1 -> 2",
    "stage2":            "Alignment & filtering (Stage 2)",
    "transition_2to3":   "Transition 2 -> 3",
    "stage3_per_family": "Per-family trees (Stage 3a)",
    "stage3_concat":     "Concatenated / species tree (Stage 3b)",
    "stage3_alerax":     "AleRax reconciliation (Stage 3c)",
}

STATUS_ICON = {
    "complete": "O",
    "skipped":  "-",
    "failed":   "X",
    "pending":  ".",
}


# ── helpers ────────────────────────────────────────────────────────────────────

def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def state_path(dataset: str) -> Path:
    return STAGE1_DIR / "resources" / dataset / ".pipeline_state.json"


def load_state(dataset: str) -> dict:
    p = state_path(dataset)
    if p.exists():
        state = json.loads(p.read_text())
        # migrate: add any steps introduced after the state file was created
        for s in STEPS:
            if s not in state["steps"]:
                state["steps"][s] = "pending"
        return state
    return {
        "dataset": dataset,
        "created": now_iso(),
        "config": {},
        "steps": {s: "pending" for s in STEPS},
        "last_updated": now_iso(),
    }


def save_state(state: dict) -> None:
    state["last_updated"] = now_iso()
    p = state_path(state["dataset"])
    p.write_text(json.dumps(state, indent=2))


def find_datasets() -> list[str]:
    resources = STAGE1_DIR / "resources"
    if not resources.exists():
        return []
    return sorted(p.parent.name for p in resources.rglob(".pipeline_state.json"))


def set_nested(d: dict, dotted_key: str, value) -> None:
    """Set a value in a nested dict using a dotted key path."""
    keys = dotted_key.split(".")
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def parse_value(s: str):
    """Try to cast a string to int or float; fall back to str."""
    for cast in (int, float):
        try:
            return cast(s)
        except ValueError:
            pass
    return s


def parse_config_with_comments(path: Path) -> list[tuple[str, str, str]]:
    """Parse a YAML config file into (dotted_key, default_value, description) tuples.

    Comment lines (# ...) on the line immediately before a key become that
    key's description.  A blank line between a comment and a key resets the
    buffer, so section-header comments are naturally ignored.
    Nested keys are joined with dots (e.g. diamond.threads).
    """
    results: list[tuple[str, str, str]] = []
    parent_stack: list[tuple[int, str]] = []  # (indent_spaces, key)
    pending: list[str] = []

    for raw in path.read_text().splitlines():
        stripped = raw.strip()

        if not stripped or stripped == "---":
            pending.clear()
            continue

        if stripped.startswith("#"):
            text = stripped[1:].strip()
            if text:
                pending.append(text)
            continue

        indent = len(raw) - len(raw.lstrip())

        # Pop parent keys that are at the same or deeper indent level.
        while parent_stack and parent_stack[-1][0] >= indent:
            parent_stack.pop()

        if ":" not in stripped:
            pending.clear()
            continue

        colon = stripped.index(":")
        key = stripped[:colon].strip()
        rest = stripped[colon + 1:].strip()

        # Strip inline comments.
        if rest.startswith("#"):
            rest = ""
        elif " #" in rest:
            rest = rest[: rest.index(" #")].strip()

        dotted = ".".join(k for _, k in parent_stack)
        dotted = f"{dotted}.{key}" if dotted else key

        if rest:
            description = " ".join(pending)
            results.append((dotted, rest, description))
        else:
            parent_stack.append((indent, key))

        pending.clear()

    return results


def rename_tree_leaves(tree_path: Path, abbrev_csv: Path, out_path: Path) -> None:
    """Replace GTDB accession ID leaf names in a Newick tree with short codes.

    Reads genome2abbrev.csv (accession,short,taxa) produced by Stage 1 and
    substitutes every accession ID found in the tree string.  Longest-first
    ordering prevents partial matches between IDs that share a prefix.
    """
    mapping: dict[str, str] = {}
    with open(abbrev_csv) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("accession"):
                continue
            parts = line.split(",")
            if len(parts) >= 2:
                mapping[parts[0].strip()] = parts[1].strip()
    tree_str = tree_path.read_text()
    for acc in sorted(mapping, key=len, reverse=True):
        tree_str = tree_str.replace(acc, mapping[acc])
    out_path.write_text(tree_str)
    print(f"  [tree] Renamed {len(mapping)} leaves: {tree_path.name} -> {out_path.name}")


# ── Stage 3 tool detection ────────────────────────────────────────────────────

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base, returning a new dict."""
    result = dict(base)
    for k, v in override.items():
        if k in result and isinstance(result[k], dict) and isinstance(v, dict):
            result[k] = _deep_merge(result[k], v)
        else:
            result[k] = v
    return result


def _load_stage3_full(cfg: dict) -> dict:
    """Return Stage 3 config.yaml defaults merged with governor cfg overrides."""
    try:
        with open(STAGE3_DIR / "config" / "config.yaml") as f:
            base = yaml.safe_load(f) or {}
    except FileNotFoundError:
        base = {}
    return _deep_merge(base, cfg)


def _binary_ok(path_str) -> bool:
    """Return True if path_str is an executable file."""
    if not path_str:
        return False
    p = Path(str(path_str))
    return p.is_file() and os.access(str(p), os.X_OK)


def _try_binary(*names: str) -> str | None:
    """Return the full path of the first binary found in PATH, or None."""
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def _module_available(module_name: str) -> bool:
    """Return True if an environment module can be shown (implies it exists)."""
    if not module_name:
        return False
    try:
        r = subprocess.run(
            ["bash", "-lc", f"module show {module_name}"],
            capture_output=True, timeout=15,
        )
        return r.returncode == 0
    except Exception:
        return False


# Tools that must be reachable before Stage 3 can run.
_STAGE3_TOOLS = [
    {
        "display":  "IQ-TREE",
        "binaries": ["iqtree3", "iqtree"],
        "path_key": ("iqtree", "path"),
        "mod_key":  ("iqtree", "module_name"),
        "set_flag": "iqtree.path",
        "mod_flag": "iqtree.module_name",
    },
    {
        "display":  "PhyloBayes-MPI",
        "binaries": ["pb_mpi"],
        "path_key": ("phylobayes", "path"),
        "mod_key":  ("phylobayes", "module_name"),
        "set_flag": "phylobayes.path",
        "mod_flag": "phylobayes.module_name",
    },
    {
        "display":  "AleRax",
        "binaries": ["alerax", "AleRax"],
        "path_key": ("alerax", "path"),
        "mod_key":  ("alerax", "module_name"),
        "set_flag": "alerax.path",
        "mod_flag": "alerax.module_name",
    },
]


def configure_stage3_tools(cfg: dict) -> dict:
    """Validate, auto-detect, and interactively configure Stage 3 tool paths.

    Merges Stage 3 config.yaml defaults with governor overrides so that a
    valid default path in config.yaml is accepted without reprompting.

    For each tool the resolution order is:
      1. Configured path exists and is executable  → use it
      2. Binary found in PATH                      → record and use it
      3. Configured module name is loadable        → use it
      4. Prompt the user (or print guidance if not on a TTY)

    OpenMPI is checked similarly: mpirun in PATH, then module, then prompt.
    Compute nodes often have OpenMPI even when the login node does not, so a
    missing module check is treated as a warning rather than a hard error.

    Returns the updated cfg dict (overrides are recorded so they persist to the
    state file and are passed to Snakemake on the next invocation).
    """
    full = _load_stage3_full(cfg)
    interactive = sys.stdin.isatty()
    unresolved: list[str] = []

    print()
    print("  Checking Stage 3 tool availability...")

    for tool in _STAGE3_TOOLS:
        k1, k2 = tool["path_key"]
        m1, m2 = tool["mod_key"]
        path_val = (full.get(k1) or {}).get(k2, "")
        mod_val  = (full.get(m1) or {}).get(m2, "")

        if _binary_ok(path_val):
            print(f"    {tool['display']:<18} found at {path_val}")
            continue

        found = _try_binary(*tool["binaries"])
        if found:
            print(f"    {tool['display']:<18} auto-detected at {found}")
            set_nested(cfg, tool["set_flag"], found)
            continue

        if _module_available(mod_val):
            print(f"    {tool['display']:<18} module '{mod_val}' available")
            continue

        # Not found — prompt or record as unresolved
        print(f"    {tool['display']:<18} NOT FOUND")
        if path_val:
            print(f"      configured path '{path_val}' does not exist or is not executable")
        if mod_val:
            print(f"      module '{mod_val}' is not loadable (may be available on compute nodes)")
        print(f"      searched PATH for: {', '.join(tool['binaries'])}")

        if interactive:
            ans = input(f"      Absolute path to {tool['display']} binary (blank to set module name instead): ").strip()
            if ans:
                set_nested(cfg, tool["set_flag"], ans)
                print(f"      Saved as {tool['set_flag']}={ans}")
            else:
                mod = input(f"      Module name for {tool['display']} (blank to skip — pipeline may fail): ").strip()
                if mod:
                    set_nested(cfg, tool["mod_flag"], mod)
                    print(f"      Saved as {tool['mod_flag']}={mod}")
                else:
                    unresolved.append(tool["display"])
        else:
            print(f"      Fix with:  --set {tool['set_flag']}=/path/to/binary")
            print(f"          or:    --set {tool['mod_flag']}=module_name")
            unresolved.append(tool["display"])

    # OpenMPI — needed by PhyloBayes-MPI and AleRax
    openmpi_mod = full.get("openmpi_module_name", "")
    if shutil.which("mpirun") or shutil.which("mpiexec"):
        print(f"    {'OpenMPI':<18} mpirun found in PATH")
    elif _module_available(openmpi_mod):
        print(f"    {'OpenMPI':<18} module '{openmpi_mod}' available")
    else:
        print(f"    {'OpenMPI':<18} not found in PATH")
        if openmpi_mod:
            print(f"      module '{openmpi_mod}' is not loadable on this node")
        print(f"      Note: compute nodes often have OpenMPI even when the login node does not.")
        if interactive:
            mod = input("      OpenMPI module name (blank = assume available on compute nodes): ").strip()
            if mod:
                cfg["openmpi_module_name"] = mod
                print(f"      Saved as openmpi_module_name={mod}")
        else:
            if openmpi_mod:
                print(f"      If this is wrong, fix with:  --set openmpi_module_name=module_name")

    if unresolved:
        print()
        print(f"  Warning: {len(unresolved)} tool(s) unresolved: {', '.join(unresolved)}")
        print("  Stage 3 may fail. Re-run interactively or supply paths via --set.")

    print()
    return cfg


# ── snakemake runner ───────────────────────────────────────────────────────────

def snakemake_unlock(stage_dir: Path) -> None:
    """Remove Snakemake lock files directly.

    After a SIGKILL, `snakemake --unlock` itself can fail because it tries to
    acquire part of the lock state. A dead process can never release them on 
    its own, so deleting files in .snakemake/locks/
    """
    locks_dir = stage_dir / ".snakemake" / "locks"
    if locks_dir.exists():
        removed = list(locks_dir.iterdir())
        for f in removed:
            f.unlink()
        if removed:
            print(f"  [unlock] Cleared {len(removed)} stale lock(s) in {locks_dir}")


def run_snakemake(
    stage_dir: Path,
    executor: str,
    jobs: int,
    dataset_name: str,
    config_overrides: dict,
    targets: list[str] | None = None,
    dry_run: bool = False,
) -> bool:
    """Invoke Snakemake in stage_dir and return True on success.

    dataset_name is always injected via --config so the correct value
    overrides whatever is in the stage's config.yaml.
    config_overrides (from --set) are written to a temp configfile that
    Snakemake merges on top of the stage's own configfile.
    """
    snakemake_unlock(stage_dir)

    cmd = [
        "snakemake",
        "--keep-incomplete",   # don't wipe partial outputs on failure
        "--rerun-incomplete",  # rerun rules whose outputs are incomplete
        "--latency-wait", "60",
        "--config", f"dataset={dataset_name}",
    ]

    if executor == "local":
        cmd += ["--cores", str(jobs)]
    else:
        cmd += ["--executor", executor, "-j", str(jobs), "--profile", "profile/"]

    if dry_run:
        cmd += ["--dry-run"]

    if targets:
        cmd += targets

    tmp_path = None
    if config_overrides:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False,
            dir=stage_dir, prefix=".governor_cfg_",
        ) as f:
            yaml.dump(config_overrides, f)
            tmp_path = Path(f.name)
        cmd += ["--configfile", str(tmp_path)]

    print(f"  $ {' '.join(str(x) for x in cmd)}")
    result = subprocess.run(cmd, cwd=stage_dir)

    if tmp_path and tmp_path.exists():
        tmp_path.unlink()

    return result.returncode == 0


# ── transition scripts ─────────────────────────────────────────────────────────

def run_transition_1to2(
    dataset: str, mcl_inflation: float, clustertype: str, clustering: str = "eggnog"
) -> bool:
    result = subprocess.run(
        ["bash", "transition_1to2.sh", dataset, str(mcl_inflation), clustertype, clustering],
        cwd=SNAKEMAKE_DIR,
    )
    return result.returncode == 0


def run_transition_2to3(dataset: str, clustering: str = "eggnog") -> bool:
    result = subprocess.run(
        ["bash", "transition_2to3.sh", dataset, clustering],
        cwd=SNAKEMAKE_DIR,
    )
    return result.returncode == 0


# ── per-dataset orchestration ──────────────────────────────────────────────────

def advance_dataset(
    dataset: str,
    executor: str,
    jobs: int,
    config_overrides: dict,
    species_tree: Path | None,
    create_db: bool,
    dry_run: bool,
) -> bool:
    state = load_state(dataset)

    # Merge new overrides into persisted config (allows changing params on resume).
    cfg = state["config"]
    cfg.update(config_overrides)
    state["config"] = cfg
    save_state(state)

    # MCL inflation and clustertype are only needed internally to locate the
    # right MCL results directory during the transition; they are not exposed
    # in dataset names.
    mcl_inflation = cfg.get("mcl_inflation", 1.8)
    clustertype = cfg.get("clustertype", "Normal")
    clustering = cfg.get("clustering", "eggnog")

    print(f"\n{'='*60}")
    print(f"Dataset: {dataset}")
    print()
    for step in STEPS:
        status = state["steps"][step]
        icon = STATUS_ICON.get(status, "?")
        print(f"  {icon}  {STEP_LABELS[step]}: {status}")
    print()

    def done(step: str) -> bool:
        return state["steps"][step] in ("complete", "skipped")

    def mark(step: str, status: str) -> None:
        state["steps"][step] = status
        save_state(state)

    def fail(step: str, msg: str) -> bool:
        mark(step, "failed")
        print(f"  [FAILED] {msg}")
        print(f"  Fix the issue and re-run this script to retry from this step.")
        return False

    # ── Stage 1 ──────────────────────────────────────────────────────────────
    if not done("stage1"):
        print(f">  {STEP_LABELS['stage1']}")
        ok = run_snakemake(STAGE1_DIR, executor, jobs, dataset, cfg, dry_run=dry_run)
        if not ok:
            return fail("stage1", "Snakemake exited non-zero.")
        mark("stage1", "complete")
    else:
        print(f"O  {STEP_LABELS['stage1']}")

    # ── Rename species tree leaves: accession IDs → short codes ─────────────
    # genome2abbrev.csv is produced by Stage 1, so this runs once Stage 1 is done.
    # The renamed tree is cached in the dataset resource directory for reruns.
    effective_species_tree: Path | None = None
    if species_tree:
        abbrev_csv = STAGE1_DIR / "results" / dataset / "genome2abbrev.csv"
        renamed = state_path(dataset).parent / "species_tree_renamed.nwk"
        if abbrev_csv.exists():
            rename_tree_leaves(species_tree, abbrev_csv, renamed)
            effective_species_tree = renamed
        else:
            print("  [tree] Warning: genome2abbrev.csv not found yet; using original leaf names")
            effective_species_tree = species_tree

    # ── Stage 1 DB: optional DuckDB/Parquet from eggnog results ──────────────
    db_status = state["steps"]["stage1_db"]
    if create_db:
        if db_status != "complete":
            print(f">  {STEP_LABELS['stage1_db']}")
            ok = run_snakemake(
                STAGE1_DIR, executor, jobs, dataset, cfg,
                targets=["db_all"], dry_run=dry_run,
            )
            if not ok:
                return fail("stage1_db", "Snakemake db_all target failed.")
            mark("stage1_db", "complete")
        else:
            print(f"O  {STEP_LABELS['stage1_db']}")
    else:
        if db_status == "complete":
            print(f"O  {STEP_LABELS['stage1_db']}")
        else:
            if db_status != "skipped":
                mark("stage1_db", "skipped")
            print(f"–  {STEP_LABELS['stage1_db']} (skipped — use --create-db to enable)")

    # ── Transition 1 → 2 ─────────────────────────────────────────────────────
    if not done("transition_1to2"):
        print(f">  {STEP_LABELS['transition_1to2']}")
        ok = run_transition_1to2(dataset, mcl_inflation, clustertype, clustering)
        if not ok:
            return fail("transition_1to2", "transition_1to2.sh failed.")
        mark("transition_1to2", "complete")
    else:
        print(f"O  {STEP_LABELS['transition_1to2']}")

    # ── Stage 2 ──────────────────────────────────────────────────────────────
    if not done("stage2"):
        print(f">  {STEP_LABELS['stage2']}")
        ok = run_snakemake(STAGE2_DIR, executor, jobs, dataset, cfg, dry_run=dry_run)
        if not ok:
            return fail("stage2", "Snakemake exited non-zero.")
        mark("stage2", "complete")
    else:
        print(f"O  {STEP_LABELS['stage2']}")

    # ── Transition 2 → 3 ─────────────────────────────────────────────────────
    if not done("transition_2to3"):
        print(f">  {STEP_LABELS['transition_2to3']}")
        ok = run_transition_2to3(dataset, clustering)
        if not ok:
            return fail("transition_2to3", "transition_2to3.sh failed.")
        mark("transition_2to3", "complete")
    else:
        print(f"O  {STEP_LABELS['transition_2to3']}")

    # ── Stage 3 tool check (runs once; skipped when all Stage 3 is done) ────────
    if not all(done(s) for s in ("stage3_per_family", "stage3_concat", "stage3_alerax")):
        cfg = configure_stage3_tools(cfg)
        state["config"] = cfg
        save_state(state)

    # ── Stage 3a: per-family trees (always required) ──────────────────────────
    if not done("stage3_per_family"):
        print(f">  {STEP_LABELS['stage3_per_family']}")
        ok = run_snakemake(
            STAGE3_DIR, executor, jobs, dataset, cfg,
            targets=["per_family_all"], dry_run=dry_run,
        )
        if not ok:
            return fail("stage3_per_family", "Snakemake exited non-zero.")
        mark("stage3_per_family", "complete")
    else:
        print(f"O  {STEP_LABELS['stage3_per_family']}")

    # ── Stage 3b: concatenated / species tree (skip if tree provided) ─────────
    if effective_species_tree:
        if not done("stage3_concat"):
            mark("stage3_concat", "skipped")
        print(f"-  {STEP_LABELS['stage3_concat']} (skipped — using provided species tree)")
        alerax_cfg = dict(cfg)
        alerax_cfg["species_tree"] = str(effective_species_tree.resolve())
    else:
        alerax_cfg = cfg
        if not done("stage3_concat"):
            print(f">  {STEP_LABELS['stage3_concat']}")
            ok = run_snakemake(
                STAGE3_DIR, executor, jobs, dataset, cfg,
                targets=["concat_all"], dry_run=dry_run,
            )
            if not ok:
                return fail("stage3_concat", "Snakemake exited non-zero.")
            mark("stage3_concat", "complete")
        else:
            print(f"O  {STEP_LABELS['stage3_concat']}")

    # ── Stage 3c: AleRax reconciliation ──────────────────────────────────────
    if not done("stage3_alerax"):
        print(f">  {STEP_LABELS['stage3_alerax']}")
        ok = run_snakemake(
            STAGE3_DIR, executor, jobs, dataset, alerax_cfg,
            targets=["alerax_all"], dry_run=dry_run,
        )
        if not ok:
            return fail("stage3_alerax", "Snakemake exited non-zero.")
        mark("stage3_alerax", "complete")
    else:
        print(f"O  {STEP_LABELS['stage3_alerax']}")

    print(f"\nO  Dataset '{dataset}' fully complete.")
    return True


# ── commands ───────────────────────────────────────────────────────────────────

def cmd_parameters(_args) -> None:
    """Print all configurable parameters parsed from each stage's config.yaml."""
    configs = [
        ("Stage 1 — Clustering", STAGE1_DIR / "config" / "config.yaml"),
        ("Stage 2 — Alignment & filtering", STAGE2_DIR / "config" / "config.yaml"),
        ("Stage 3 — Phylogenetic inference", STAGE3_DIR / "config" / "config.yaml"),
    ]
    # Keys that are pipeline-internal and not meaningful to override via --set.
    skip = {"dataset", "accession_ids_file", "manually_add_taxa_file"}

    for label, path in configs:
        if not path.exists():
            continue
        params = [
            (k, v, d)
            for k, v, d in parse_config_with_comments(path)
            if k not in skip
        ]
        if not params:
            continue
        print(f"\n{label}")
        print("-" * len(label))
        for key, default, desc in params:
            desc_str = f"  {desc}" if desc else ""
            print(f"  {key:<45} (default: {default}){desc_str}")

    print()
    print("Pass overrides as:  ./run_pipeline.py run --set key=value [key=value ...]")
    print("Nested keys use dot notation:  --set diamond.threads=64 phylobayes.tasks=32")


def cmd_status(_args) -> None:
    datasets = find_datasets()
    if not datasets:
        print("No datasets found. Start one with:  run_pipeline.py run --dataset NAME --ids FILE")
        return
    for ds in datasets:
        state = load_state(ds)
        steps = state["steps"]
        n_done = sum(1 for s in steps.values() if s in ("complete", "skipped"))
        print(f"\n{ds}  [{n_done}/{len(steps)} steps done]  last updated: {state['last_updated']}")
        for step in STEPS:
            status = steps[step]
            icon = STATUS_ICON.get(status, "?")
            print(f"  {icon}  {STEP_LABELS[step]}: {status}")


def cmd_run(args) -> None:
    # Build config overrides from --set key=value args
    config_overrides: dict = {}
    for kv in (args.set or []):
        if "=" not in kv:
            print(f"Error: --set value must be key=value, got: {kv!r}", file=sys.stderr)
            sys.exit(1)
        k, _, v = kv.partition("=")
        set_nested(config_overrides, k.strip(), parse_value(v.strip()))

    if args.shared_data_dir:
        config_overrides["shared_data_dir"] = args.shared_data_dir

    if args.clustering:
        config_overrides["clustering"] = args.clustering

    if args.abbrev_map:
        abbrev_map_path = Path(args.abbrev_map)
        if not abbrev_map_path.exists():
            print(f"Error: abbreviation map file not found: {abbrev_map_path}", file=sys.stderr)
            sys.exit(1)
        config_overrides["custom_abbrev_map"] = str(abbrev_map_path.resolve())

    species_tree = None
    if args.species_tree:
        species_tree = Path(args.species_tree)
        if not species_tree.exists():
            print(f"Error: species tree file not found: {species_tree}", file=sys.stderr)
            sys.exit(1)

    # Resolve which datasets to process
    if args.dataset:
        datasets = [args.dataset]
        dd = STAGE1_DIR / "resources" / args.dataset
        dd.mkdir(parents=True, exist_ok=True)
        if args.ids:
            shutil.copy(args.ids, dd / "list_of_accession_ids.txt")
        elif not (dd / "list_of_accession_ids.txt").exists():
            print(
                f"Error: no accession ID list found for dataset '{args.dataset}'.\n"
                f"Provide --ids FILE to supply one.",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        if args.ids:
            print("Error: --ids requires --dataset to be specified.", file=sys.stderr)
            sys.exit(1)
        datasets = find_datasets()
        if not datasets:
            print(
                "No datasets found. Start one with:\n"
                "  run_pipeline.py run --dataset NAME --ids FILE",
                file=sys.stderr,
            )
            sys.exit(1)

    any_failed = False
    for ds in datasets:
        ok = advance_dataset(
            dataset=ds,
            executor=args.executor,
            jobs=args.jobs,
            config_overrides=dict(config_overrides),
            species_tree=species_tree,
            create_db=args.create_db,
            dry_run=args.dry_run,
        )
        if not ok:
            any_failed = True

    sys.exit(1 if any_failed else 0)


# ── CLI ────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Phylogenomics pipeline governor.\n"
            "Invoke at any time to check and advance all datasets from where they left off."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # status
    sub.add_parser("status", help="Show current state of all known datasets.")

    # parameters
    sub.add_parser(
        "parameters",
        help="List all configurable parameters with their defaults and descriptions.",
    )

    # run
    r = sub.add_parser(
        "run",
        help="Start or resume the pipeline for one or all datasets.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "If --dataset is omitted, all known datasets are advanced.\n"
            "Re-running after a failure resumes from the failed step.\n\n"
            "Examples:\n"
            "  Start a new dataset:\n"
            "    run_pipeline.py run --dataset my_run --ids genomes.txt\n\n"
            "  Resume all datasets after a SLURM failure:\n"
            "    run_pipeline.py run\n\n"
            "  Use a pre-existing species tree:\n"
            "    run_pipeline.py run --dataset my_run --species-tree tree.nwk\n\n"
            "  Use bulk storage for large files (Deigo /bucket layout):\n"
            "    run_pipeline.py run --dataset my_run --ids genomes.txt \\\n"
            "        --shared-data-dir /bucket/user/project\n\n"
            "  Override config params:\n"
            "    run_pipeline.py run --dataset my_run --set eggnog.threads=32 iqtree.threads=64\n"
        ),
    )
    r.add_argument(
        "--dataset", metavar="NAME",
        help="Dataset name. Omit to advance all known datasets.",
    )
    r.add_argument(
        "--ids", metavar="FILE",
        help="File of GTDB accession IDs, one per line (required for new datasets).",
    )
    r.add_argument(
        "--executor", default="slurm", metavar="EXECUTOR",
        help=(
            "Snakemake executor plugin (e.g. slurm, local, etc). "
            "Default: slurm."
        ),
    )
    r.add_argument(
        "--jobs", type=int, default=20, metavar="N",
        help="Maximum parallel jobs / cores passed to Snakemake. Default: 20.",
    )
    r.add_argument(
        "--species-tree", metavar="FILE",
        help=(
            "Path to a pre-existing species tree in Newick format. "
            "When provided, the concatenated inference step (Stage 3b) is skipped "
            "and this tree is passed directly to AleRax."
            "Leaf names of spcies tree should match GTDB genome accession IDs."
        ),
    )
    r.add_argument(
        "--set", nargs="*", metavar="key=value",
        help=(
            "Override any pipeline config parameter using dotted key paths. "
            "Values are auto-cast to int/float where possible. "
            "Example: --set eggnog.threads=32 phylobayes.tasks=64. "
            "Advanced: use mcl_inflation=1.4 or clustertype=Normal to select "
            "which MCL result directory feeds into Stage 2 (default: 1.8 / Normal)."
        ),
    )
    r.add_argument(
        "--clustering", choices=["eggnog", "mcl"], default=None,
        help=(
            "Clustering method for gene family assignment. "
            "'eggnog' (default): only COG-based clustering via EggNOG-mapper; "
            "skips DIAMOND all-vs-all and MCL entirely, saving substantial compute. "
            "'mcl': also runs MCL with all configured inflation factors in parallel "
            "alongside EggNOG, useful for comparison or when EggNOG annotation is "
            "unavailable. Stored in the dataset state once set."
        ),
    )
    r.add_argument(
        "--abbrev-map", metavar="FILE",
        help=(
            "CSV file mapping GTDB accession IDs to custom leaf-node names "
            "(columns: accession,short). When provided, this replaces the "
            "auto-generated abbreviations throughout the pipeline. "
            "If --species-tree is also given, its leaves are expected to use "
            "GTDB accession IDs and will be renamed to match."
        ),
    )
    r.add_argument(
        "--shared-data-dir", metavar="PATH",
        help=(
            "Path to bulk storage for large shared files (GTDB archives, extracted "
            "genome directories, taxonomy/metadata TSVs, eggnog DB, and protein "
            "parquets from the DB step). It's useful if the machine has a small "
            "execution storage (i.e. /scratch), but a large backstorage that is "
            "slow to write, but fast enough to read. "
            "Stored in the dataset state -- only needs to be specified once. "
            "When unset, all files stay under the project resources/ directory."
        ),
    )
    r.add_argument(
        "--create-db", action="store_true",
        help=(
            "After Stage 1, build a DuckDB + Parquet database from the eggnog-mapper "
            "results. Produces results/{dataset}/db/dataset.duckdb "
            "and two Parquet files (genome.parquet, protein.parquet) that can be queried "
            "with SQL via DuckDB. Requires duckdb and pyarrow in the environment."
        ),
    )
    r.add_argument(
        "--dry-run", action="store_true",
        help="Pass --dry-run to Snakemake: show what would be done without running.",
    )

    args = parser.parse_args()
    if args.command == "status":
        cmd_status(args)
    elif args.command == "parameters":
        cmd_parameters(args)
    else:
        cmd_run(args)


if __name__ == "__main__":
    main()
