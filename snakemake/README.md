# Snakemake pipeline

## Clone pipeline repository

Clone this repository (in case of HPC application: into a storage space that is writable by the compute nodes):
```
cd /to-your-work-directory

git clone https://github.com/oist/subarashii-annotation-pipeline.git
cd subarashii-annotation-pipeline/
git checkout -b snakemake origin/snakemake
cd snakemake
```

## Installing Snakemake

Installing Minforge (formerly known as Mambaforge):
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# during the instlattaion below, choose an installation directory writable in 
# an HPC setting and with enough storage space to install future dependencies
bash Miniforge3-$(uname)-$(uname -m).sh
# to make the changes take effect, logout and login or execute:
source ~/.bashrc
```

Create environment with Snakemake:
```
# manual creation, if you choose this you have to install packages manually later
#conda create -c conda-forge -c bioconda -n snakemake_env snakemake
# see env.yml for packages, installation command:
#mamba install -c conda-forge -c bioconda <packagename>

# creating and fetching necessary dependencies automatically:
mamba create -f env.yml
```

Activate environment:
```
mamba activate snakemake_env
```

## Preparing your dataset

Choose a short and descriptive name for your dataset without any special characters and create a directory for it:
```
mkdir -p 1_cluster/resources/my_awesome_dataset
```

Also set this name in the `1_cluster/config/config.yaml`, `2_concatenate_and_filter/config/config.yaml` and `3_inference/config/config.yaml`:
```
dataset: "my_awesome_dataset"
```

Create a file named `list_of_accession_ids.txt` inside the directory `1_cluster/resources/my_awesome_dataset` that has one accession id per line and then a newline character, e.g.:
```
RS_GCF_000020505.1
GB_GCA_004322215.1
RS_GCF_000419585.1
GB_GCA_001803045.1
```

Decide which GTDB version you want to use and modify `1_cluster/config/config.yaml` accordingly, e.g. for GTDB v89.0, you would set:
```
gtdb_version: "89"
gtdv_subversion: "0"
```

Some accession IDs may not be found in GTDB, but they are in NCBI (in which case the pipeline will download it from NCBI), but then the GTDB's taxonomy file won't contain metadata about these genomes. To add genome metadata manually (in the same format as the GTDB taxonomy file), you may use the `manually_add_taxa_file` variable in the `1_cluster/config/config.yaml`.

## Personalizing settings

Some use-case and environment-related settings are collected in the `config/config.yml` file.

Edit it with focus on the variables `evalue_cut`, `cov_cut`, `runtime`, etc. to fit to your case.

## How to run snakemake in general (conceptually, don't do it yet)

If running on remote server, start snakemake from inside a virtual temrinal like `screen` or `tmux`.

With slurm:
```
snakemake --executor slurm -j 20 --profile profile/ --slurm-keep-successful-logs
```

Locally:
```
snakemake --cores 1
```

## Running the pipeline

The pipeline is split into 3 subpipeline, due to the way Snakemake works (it has to know the names of the input files, so we cannot have rules doing something with e.g. the gene families, when we don't even know yet, how many gene families we have and what are their names):
- `1_cluster`
- `2_concatenate_and_filter`
- `3_inference`

Between the subpipleines, axularily scripts `transfer_1to2.sh` and `transfer_2to3.sh` are helping to copy the results of one subpipline into the resources directory of the next one.

So first enter the directory `1_cluster`. Let's assume we have our projectname `awesome`. Create `1_cluster/resources/awesome` directory and a file name `list_of_accession_ids.txt` inside this directory. The `txt` file should contain a GTDB accession ID each line.

Next we have to customize the Snakemake settings as described above, i.e. setting the GTDB release and other parameters.

Once done, we can execute Snakemake as discussed above (either with slurm or locally, settings the cores and job numbers appropriately).
