
# Run the pipeline

## Clone the repository

```bash
git clone https://github.com/theislab/spapros-smk.git
```

## Set up snakemake environment

(similar to [snakemake install](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html))

First install mamba in your base environment in case you don't have it already:

```bash 
conda install -n base -c conda-forge mamba

```

Then create a new environment with snakemake:

```bash
conda activate base
mamba create -c conda-forge -c bioconda -n spapros_smk snakemake
```

## Set experiment configuration and test as dry-run

To set the config edit the file `spapros-smk/config/config.yaml`. (details see below)

```bash
conda activate spapros_smk
cd spapros-smk

snakemake -n c1??
```

## Run on current machine

E.g. using 8 cores

```bash
snakemake --cores 8 --use-conda --use-singularity
```

## Run on cluster

The following describes how to run the pipeline on a cluster (using the [slurm](https://slurm.schedmd.com/documentation.html) workload manager). To run it on your cluster you might need/want to adjust the configuration file `spapros-smk/config/cluster/config.yaml`.

In the `spapros-smk` directory create a file, e.g. `snakemake.sbatch` with the following content (adjust partition according to your cluster):

```bash
#!/bin/bash
  
#SBATCH -o logs/snakemake_output_%j.txt
#SBATCH -e logs/snakemake_output_%j.txt
#SBATCH -J sp_snakemake
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 3-00:00:00
#SBATCH --nice=0
#SBATCH -p cpu_p

conda activate spapros_smk

export SINGULARITY_TMPDIR=<some/tmp/dir/that/exists/and/has/space>
export APPTAINER_TMPDIR=<some/tmp/dir/that/exists/and/has/space>

snakemake --profile config/cluster/

```

To run the full pipeline with jobs distributed over the cluster submit the main snakemake job (it continues after disconnecting from the cluster):

```bash
sbatch snakemake.sbatch
```


(Alternatively) It's also possible to just run the snakemake command with the cluster profile without submitting an extra job:
```bash
snakemake --profile config/cluster/
```

## Important notes to not introduce errors

- Do not delete the files `data_parameters.csv`, `evaluation_overview.csv`, `selection_overview.csv` and `selection_parameters.csv` in case you change the config.yaml for additional runs. They will be extended automatically, otherwise the new ids in the new csvs might not match to the old ones and then incorrect datasets etc. are loaded.
- Extending the `config.yaml` with new runs is fine. Also old runs can be deleted from the config.yaml. If they were already calculated, the old results aren't deleted. 
- If this is confusing, you can also delete everything in the data_tmp and results folders and restart the runs. This will take longer, but you can be sure that the results are correct.


## Configuration file
TODO: add description of config file