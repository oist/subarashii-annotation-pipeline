#!/usr/bin/env python3
"""
Phylogenomics pipeline governor.

Run at any time — it checks where each dataset is and continues from there.
State is persisted between invocations, so re-running after a failure resumes
from the failed step without redoing completed work.

Usage examples:
  # Start a new dataset
  ./run_pipeline.py run --dataset my_run --ids genomes.txt

  # Resume all known datasets (e.g. after a SLURM failure)
  ./run_pipeline.py run

  # Check status without running anything
  ./run_pipeline.py status

  # Use a pre-existing species tree (skips concatenated inference step)
  ./run_pipeline.py run --dataset my_run --species-tree /path/to/tree.nwk

  # Override config parameters
  ./run_pipeline.py run --dataset my_run --executor local --jobs 8 --set eggnog.threads=16
"""

import argparse
import json
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
    "transition_1to2",
    "stage2",
    "transition_2to3",
    "stage3_per_family",
    "stage3_concat",
    "stage3_alerax",
]

STEP_LABELS = {
    "stage1":           "Clustering (Stage 1)",
    "transition_1to2":  "Transition 1 → 2",
    "stage2":           "Alignment & filtering (Stage 2)",
    "transition_2to3":  "Transition 2 → 3",
    "stage3_per_family": "Per-family trees (Stage 3a)",
    "stage3_concat":    "Concatenated / species tree (Stage 3b)",
    "stage3_alerax":    "AleRax reconciliation (Stage 3c)",
}

STATUS_ICON = {
    "complete": "✓",
    "skipped":  "–",
    "failed":   "✗",
    "pending":  "·",
}


# ── helpers ────────────────────────────────────────────────────────────────────

def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def inflated_name(dataset: str, inflation: float, clustertype: str = "Normal") -> str:
    """Compute the dataset name used by Stages 2 and 3.

    Mirrors what transition_1to2.sh produces:
      {dataset}_I{inflation_without_dot}_{clustertype}
    e.g. my_data + 1.8 + Normal → my_data_I18_Normal
    """
    inf_str = str(inflation).replace(".", "")
    return f"{dataset}_I{inf_str}_{clustertype}"


def state_path(dataset: str) -> Path:
    return STAGE1_DIR / "resources" / dataset / ".pipeline_state.json"


def load_state(dataset: str) -> dict:
    p = state_path(dataset)
    if p.exists():
        return json.loads(p.read_text())
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


# ── snakemake runner ───────────────────────────────────────────────────────────

def snakemake_unlock(stage_dir: Path) -> None:
    subprocess.run(["snakemake", "--unlock"], cwd=stage_dir, capture_output=True)


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

def run_transition_1to2(dataset: str, inflation: float, clustertype: str) -> bool:
    result = subprocess.run(
        ["bash", "transition_1to2.sh", dataset, str(inflation), clustertype],
        cwd=SNAKEMAKE_DIR,
    )
    return result.returncode == 0


def run_transition_2to3(inf_dataset: str) -> bool:
    result = subprocess.run(
        ["bash", "transition_2to3.sh", inf_dataset],
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
    dry_run: bool,
) -> bool:
    state = load_state(dataset)

    # Merge new overrides into persisted config (allows changing params on resume).
    cfg = state["config"]
    cfg.update(config_overrides)
    state["config"] = cfg
    save_state(state)

    inflation = cfg.get("mcl_inflation", 1.8)
    clustertype = cfg.get("clustertype", "Normal")
    inf_dataset = inflated_name(dataset, inflation, clustertype)

    print(f"\n{'='*60}")
    print(f"Dataset : {dataset}")
    print(f"Stage 2/3 name : {inf_dataset}")
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
        print(f"▶  {STEP_LABELS['stage1']}")
        ok = run_snakemake(STAGE1_DIR, executor, jobs, dataset, cfg, dry_run=dry_run)
        if not ok:
            return fail("stage1", "Snakemake exited non-zero.")
        mark("stage1", "complete")
    else:
        print(f"✓  {STEP_LABELS['stage1']}")

    # ── Transition 1 → 2 ─────────────────────────────────────────────────────
    if not done("transition_1to2"):
        print(f"▶  {STEP_LABELS['transition_1to2']}")
        ok = run_transition_1to2(dataset, inflation, clustertype)
        if not ok:
            return fail("transition_1to2", "transition_1to2.sh failed.")
        mark("transition_1to2", "complete")
    else:
        print(f"✓  {STEP_LABELS['transition_1to2']}")

    # ── Stage 2 ──────────────────────────────────────────────────────────────
    if not done("stage2"):
        print(f"▶  {STEP_LABELS['stage2']}")
        ok = run_snakemake(STAGE2_DIR, executor, jobs, inf_dataset, cfg, dry_run=dry_run)
        if not ok:
            return fail("stage2", "Snakemake exited non-zero.")
        mark("stage2", "complete")
    else:
        print(f"✓  {STEP_LABELS['stage2']}")

    # ── Transition 2 → 3 ─────────────────────────────────────────────────────
    if not done("transition_2to3"):
        print(f"▶  {STEP_LABELS['transition_2to3']}")
        ok = run_transition_2to3(inf_dataset)
        if not ok:
            return fail("transition_2to3", "transition_2to3.sh failed.")
        mark("transition_2to3", "complete")
    else:
        print(f"✓  {STEP_LABELS['transition_2to3']}")

    # ── Stage 3a: per-family trees (always required) ──────────────────────────
    if not done("stage3_per_family"):
        print(f"▶  {STEP_LABELS['stage3_per_family']}")
        ok = run_snakemake(
            STAGE3_DIR, executor, jobs, inf_dataset, cfg,
            targets=["per_family_all"], dry_run=dry_run,
        )
        if not ok:
            return fail("stage3_per_family", "Snakemake exited non-zero.")
        mark("stage3_per_family", "complete")
    else:
        print(f"✓  {STEP_LABELS['stage3_per_family']}")

    # ── Stage 3b: concatenated / species tree (skip if tree provided) ─────────
    if species_tree:
        if not done("stage3_concat"):
            mark("stage3_concat", "skipped")
        print(f"–  {STEP_LABELS['stage3_concat']} (skipped — using provided species tree)")
        alerax_cfg = dict(cfg)
        alerax_cfg["species_tree"] = str(species_tree.resolve())
    else:
        alerax_cfg = cfg
        if not done("stage3_concat"):
            print(f"▶  {STEP_LABELS['stage3_concat']}")
            ok = run_snakemake(
                STAGE3_DIR, executor, jobs, inf_dataset, cfg,
                targets=["concat_all"], dry_run=dry_run,
            )
            if not ok:
                return fail("stage3_concat", "Snakemake exited non-zero.")
            mark("stage3_concat", "complete")
        else:
            print(f"✓  {STEP_LABELS['stage3_concat']}")

    # ── Stage 3c: AleRax reconciliation ──────────────────────────────────────
    if not done("stage3_alerax"):
        print(f"▶  {STEP_LABELS['stage3_alerax']}")
        ok = run_snakemake(
            STAGE3_DIR, executor, jobs, inf_dataset, alerax_cfg,
            targets=["alerax_all"], dry_run=dry_run,
        )
        if not ok:
            return fail("stage3_alerax", "Snakemake exited non-zero.")
        mark("stage3_alerax", "complete")
    else:
        print(f"✓  {STEP_LABELS['stage3_alerax']}")

    print(f"\n✓  Dataset '{dataset}' fully complete.")
    return True


# ── commands ───────────────────────────────────────────────────────────────────

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

    if args.inflation is not None:
        config_overrides["mcl_inflation"] = args.inflation
    if args.clustertype is not None:
        config_overrides["clustertype"] = args.clustertype

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
            "Snakemake executor plugin (e.g. slurm, local, lsf, google-lifesciences). "
            "Default: slurm."
        ),
    )
    r.add_argument(
        "--jobs", type=int, default=20, metavar="N",
        help="Maximum parallel jobs / cores passed to Snakemake. Default: 20.",
    )
    r.add_argument(
        "--inflation", type=float, metavar="FLOAT",
        help=(
            "MCL inflation value to use for the Stage 1→2 transition. "
            "Must match one of the values in Stage 1's config. Default: 1.8."
        ),
    )
    r.add_argument(
        "--clustertype", metavar="STR",
        help="Cluster type used in the transition (Normal / …). Default: Normal.",
    )
    r.add_argument(
        "--species-tree", metavar="FILE",
        help=(
            "Path to a pre-existing species tree in Newick format. "
            "When provided, the concatenated inference step (Stage 3b) is skipped "
            "and this tree is passed directly to AleRax."
        ),
    )
    r.add_argument(
        "--set", nargs="*", metavar="key=value",
        help=(
            "Override any pipeline config parameter using dotted key paths. "
            "Values are auto-cast to int/float where possible. "
            "Example: --set eggnog.threads=32 phylobayes.tasks=64"
        ),
    )
    r.add_argument(
        "--dry-run", action="store_true",
        help="Pass --dry-run to Snakemake: show what would be done without running.",
    )

    args = parser.parse_args()
    if args.command == "status":
        cmd_status(args)
    else:
        cmd_run(args)


if __name__ == "__main__":
    main()
