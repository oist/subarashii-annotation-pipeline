# Snakemake pipeline

## Clone pipeline repository

Clone this repository (in case of HPC application: into a storage space that is writable by the compute nodes):
```
cd /workpath
git clone git@github.com:oist/subarashii-annotation-pipeline.git
cd subarashii-annotation-pipeline/snakemake
```

## Installing Snakemake

Installing Minforge (formerly known as Mambaforge):
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
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
conda activate snakemake_env
```

## Personalizing settings

Some use-case and environment-related settings are collected in the `config/config.yml` file.

Edit it with focus on the variables `` and ``, to fit your case.

## Running pipeline

With slurm:
```
snakemake --executor slurm -j 100 --profile profile/
```

Locally:
```
snakemake --cores 1
```

