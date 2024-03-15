# pipeline_rnaseq_rte

- Create a project directory
- In the same parent directory as project, clone this pipeline_rnaseq_rte directory
- Move the contents of the move_to_project_dir_and_edit directory to your project directory.
- Make all of the environments in the envs/ directory. E.g.:
   ```
   mamba env create --file envs/rseqc.yaml
   ```
- Modify the contents of conf/private/sample_table. Make sure sample names do not start with numbers; add an X in front if they do.
- Modify the contents of conf/private/configPrivate.yaml; make sure the samples, contrast, library_type are modified properly
- Make sure the contents of conf/shared/configShared.yaml are all paths you have access to and don't give you permissions errors
  ```
  cat {path}
  ```
- Create, in your project directory, the rawdata directory, and move your fastqs there. Make sure the naming is consistent with the naming scheme set forth in the conf/private/project_config.yaml i.e.:

```
source1: "rawdata/{sample_name}_R1.fastq.gz" source2: "rawdata/{sample_name}_R2.fastq.gz"
```

- Don't worry about peptable.csv, this is automatically updated each time you call snakemake.
- In the project_dir snakefile edit the following path so that it reflects where the pipeline snakefile lives:

```
    module mainworkflow:
            snakefile: "../pipeline_rnaseq_rte/snakefile"
```

- Adjust all resource / threads parameters according to the maximum you have availible on your system; i.e threads: 32 -> threads: 10
- Create a convenience alias for calling snakemake and sending jobs to slurm: add this to your ~/.bashrc or ~/.zshrc:

```
alias ss="snakemake \
    --software-deployment-method conda \
    --conda-frontend mamba \
    --rerun-triggers mtime \
    --keep-going \
    --rerun-incomplete  \
    --latency-wait 30 \
    --executor slurm \
    --default-resources \
        mem_mb=30000 \
        disk_mb=200000 \
        slurm_account=default \
        slurm_partition=batch \
        threads=1 \
        runtime=300 \
    --jobs 30"
```

- Install snakemake executor plugins in your snamemake conda environment, e.g.
    ```
    conda activate your_snakemake_env
    mamba install snakemake-executor-plugin-slurm
    ```
