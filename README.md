# pipeline_rnaseq_rte

Copy the example_conf directory to your project directory, and rename it to "conf".
Modify the contents of conf/private/sample_table. Make sure sample names do not start with numbers; add an X in front if they do.
Modify the contents of conf/private/configPrivate.yaml:
  make sure the samples, contrast, library_type
Create the rawdata directory, and move your fastqs there. Make sure the naming is consistent with the naming scheme set forth in the conf/private/project_config.yaml i.e.:       source1: "rawdata/{sample_name}_R1.fastq.gz" source2: "rawdata/{sample_name}_R2.fastq.gz"
Don't worry about peptable.csv, this is automatically updated each time you call snakemake.

In the project_dir snakefile edit the following path so that it reflects where the pipeline snakefile lives:
```
    module mainworkflow:
            snakefile: "../pipeline_rnaseq/snakefile"
```

In the pipeline snakefile, replace all conda environment names with a full path to Maxfield's conda environment
    i.e conda: "deseq" -> conda: "/users/mkelsey/anaconda/deseq"

Also, adjust all resource / threads paraters according to the maximum you have availible on your system.
    i.e threads: 32 -> threads: 10


Create a convenience alias for calling snakemake and sending jobs to slurm:
Add this to your ~/.bashrc or ~/.zshrc:
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
