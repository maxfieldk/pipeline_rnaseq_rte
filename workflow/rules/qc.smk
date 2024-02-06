rule mycoplasmaCheck:
    input:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
    output:
        sam = "qc/mycoplasma/mycoplasma{sample}.sam"
    threads: 12
    resources:
        mem_mb  = 30000
    conda:
        "omics"
    shell:
        """
mkdir -p $(dirname {output.sam})
bowtie2 --threads 10 -x /users/mkelsey/data/ref/genomes/mycoplasma/mycoplasma_index -1 {input.r1} -S {output.sam}
samtools stats {output.sam} > {output.sam}.stats.txt
        """


checkpoint inferLibraryType:
    input:
        bam = expand("outs/{sample}/star_output/{sample}.sorted.primary.bam", sample = samples[0])
    output:
        librarytype = "qc/library_type.txt"
    params:
        gtf = config['annotation_genes_bed12'],
    resources:
        mem_mb  = 30000
    conda:
        "rseqc"
    shell:
        """
infer_experiment.py -r {params.gtf} -i {input.bam} > {output.librarytype}
        """

rule multiqc:
    input:
        bams = expand("outs/{sample}/star_output/{sample}.sorted.bam.stats.txt", sample = samples),
        deseq = expand("results/agg/deseq2/{counttype}outfile.txt", counttype = counttypes)
    output:
        "qc/multiqc/multiqc_report.html"
    conda:
        "qc"
    shell:
        """
mkdir -p qc/multiqc
multiqc -f -o qc/multiqc --export ./
        """



