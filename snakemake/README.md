# Snakemake pipeline

## Installing Snakemake

Installing Mumbaforge:
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
source ~/.bashrc
```

Create environment wiht Snakemake:
```
conda create -c conda-forge -c bioconda -n snakemake_env snakemake
conda activate snakemake_env
```

## Running pipeline

```
snakemake --cores 1
```
