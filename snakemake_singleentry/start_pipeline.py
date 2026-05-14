#!/usr/bin/env python3
"""
Subarashii Annotation Pipeline (SAP) starter script

The SAP allows one to perform all the usual steps of a phylogenetic
inference workflow in an unattended manner.

The SAP will customize each step to fit your data size and type, so
that you just have to come back for the results and can invest your
time into research, instead of writing scripts.

Authors:
    - Lenard L. Szantho <lenardszantho@gmail.com>
    - Adrian A. Davin <aa@gmail.com>

Version:
    - v0.1 initial release (2026-01-14)

"""
import os
import sys
import argparse
from pathlib import Path
import subprocess

import snakemake
import yaml

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("project_name", help="Unique identifier of the project during the annotation process")
    parser.add_argument("gtdb_accession_ids", help="List of GTDB accession IDs or a path to a file containing the list of GTDB accession IDs to be processed by the Subarashii Annotation Pipeline", nargs="+")

    parser.add_argument("--gtdb-version", help="Specify the GTDB release to be used, by default r226 (2025) is used", default="226")
    parser.add_argument("--slurm-partition", help="Set the name of the Slurm partition Subarashii Annotation Pipeline should use for executing the tasks", default="compute")
    parser.add_argument("--slurm-maxcpus", help="Set the maximum number of cores a typical compute node has on the Slurm cluster", default=64)

    args = parser.parse_args()

    workflow_dir = "1_cluster"
    
    # read accession IDs
    gtdb_acc_ids = []
    if len(args.gtdb_accession_ids) > 1:
        # it's a list of accession IDs, not a file
        for ids in args.gtdb_accession_ids:
            gtdb_acc_ids.append(ids.strip())
    else:
        # doesn't make sense to run pipeline for a single accession ID
        if os.path.exists(args.gtdb_accession_ids[0]):
            with open(args.gtdb_accession_ids[0], "r") as inputfh:
                for rows in inputfh:
                    gtdb_acc_ids.append(rows.strip())
        else:
            print(f"File of accession IDs {args.gtdb_accession_ids} does not exist.")
            sys.exit(102)


    number_of_accession_ids = len(gtdb_acc_ids)

    Path(f"{workflow_dir}/resources/{args.project_name}").mkdir(parents=True, exist_ok=True)
    with open(f"{workflow_dir}/resources/{args.project_name}/list_of_accession_ids.txt", "w") as outputfh:
        for ids in gtdb_acc_ids:
            outputfh.write(f"{ids}\n")

    # prepare profile file

    profile_description = {
        'default-resources': {
            'slurm_partition': "compute",
            'mem_mp_per_cpu': 2000,
            'runtime': 120,
            'tasks': 1
            },
        'set-resources': {
            'get_all_proteomes': {
                'runtime': 1440
                },
            'get_all_genomes': {
                'runtime': 1440
                },
            'select_proteomes': {
                'runtime': 1440
                },
            'create_genome2taxa': {
                'runtime': 600
                },
            'make_diamond_db': {
                'runtime': 600
                },
            'diamond_all_vs_all': {
                'runtime': 5760,
                'mem_mb_per_cpu': 10000
                },
            'rename_proteins': {
                'runtime': 600
                },
            'eggnog_mapper': {
                'runtime': 5760
                }
           }
    }

    Path(f"{workflow_dir}/profile").mkdir(parents=True, exist_ok=True)
    with open(f"{workflow_dir}/profile/config.yaml", 'w') as outputfh:
        yaml.dump(profile_description, outputfh, default_flow_style=True)

    # prepare config file

    config_description = {
        'dataset': args.project_name,
        'accession_ids_file': "list_of_accession_ids.txt",
        'gtdb_version': "226",
        'gtdb_subversion': "0",
        'gtdb_base_url': "https://data.ace.uq.edu.au/public/gtdb/data/releases/",
        'manually_add_taxa_file': "manually_add_taxa.txt",
        'diamond': {
            'threads': 32,
            'evalue_cut': 1e-5,
            'cov_cut': 50 # percentage
            },
        'mcl': {
            'inflation': [1.4, 1.8, 2.0],
            'memory': 100000,
            'runtime': 1440
            },
        'eggnog': {
            'threads': 60
            },
        'mafft': {
            'threads': 4
            }
    }

    Path(f"{workflow_dir}/config").mkdir(parents=True, exist_ok=True)
    with open(f"{workflow_dir}/config/config.yaml", 'w') as outputfh:
        yaml.dump(config_description, outputfh, default_flow_style=True)    

    current_dir = os.path.dirname(os.path.realpath(__file__)) 
    snakemake_args = [
            "snakemake",
            "--executor", "slurm",
            "-j" "20",
            #"--unlock",
            "--slurm-keep-successful-logs",
            f"--snakefile={workflow_dir}/Snakefile",
            f"--directory={workflow_dir}",
            f"--profile={workflow_dir}/profile/",
            ]

    #sm_exit = snakemake.main(snakemake_args)
    sm_exit = subprocess.run(snakemake_args, text=True)
    # env=os.environ.copy(), stdout=sys.stdout, stderr=sys.stderr)

    #print(sm_exit.stderr.strip())
    #print(sm_exit.stdout.strip())
    if sm_exit != 0:
        print("Some error happened")
     

if __name__ == "__main__":
    main()
