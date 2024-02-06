rule prefetch:
    output:
        sra = temp("rawdata/{sample}/{sample}.sra")
    threads: 4
    conda:
        "omics"
    shell: "prefetch {wildcards.sample} --output-directory rawdata"

rule fastqdump:
    input: "rawdata/{sample}/{sample}.sra"
    threads: 4
    params:
        outdir = "rawdata/{sample}"
    log: "logs/{sample}/fastqdump.log"
    output:
        r1 = temp("rawdata/{sample}_1.fastq"),
        r2 = temp("rawdata/{sample}_2.fastq")
    conda:
        "omics"
    shell: "fastq-dump --split-files --outdir {params.outdir} {input} 2> {log}"

rule fastp:
    input:
        r1=lambda wildcards: peptable.loc[peptable["sample_name"] == wildcards.sample, "R1"].iloc[0],
        # r2=lambda wildcards: peptable.loc[peptable["sample_name"] == wildcards.sample, "R2"].iloc[0]
    priority: 100
    threads: 6
    log: "logs/{sample}/fastp.log"
    output:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        # r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz",
        json = "outs/{sample}/trimmedReads/fastp.json",
        html = "outs/{sample}/trimmedReads/fastp.html"
    conda:
        "qc"
    shell:
        """
fastp -i {input.r1} -o {output.r1} --json {output.json} --html {output.html} --thread {threads}
        """
