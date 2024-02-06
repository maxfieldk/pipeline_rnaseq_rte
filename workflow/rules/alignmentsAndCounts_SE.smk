
rule alignSTAR:
    input:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        # r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz"
    params:
        index = config["starindex"],
        outdir = "outs/{sample}/star_output/"
    output:
        bam = temp("outs/{sample}/star_output/Aligned.out.sam")
    threads: 8
    resources:
        mem_mb  = 60000
    conda:
        "star"
    shell:
        """
STAR --genomeDir {params.index} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {params.outdir} --runThreadN {threads} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100
        """



rule featurecounts_genes:
    input:
        sortedSTARbams = expand("outs/{sample}/star_output/{sample}.sorted.primary.bam", sample = samples),
        libtype = "qc/library_type.txt"
    output:
        countsmessy = "outs/agg/featurecounts_genes/counts_messy.txt",
        counts = "outs/agg/featurecounts_genes/counts.txt",
        countsstrandnonspecificmessy = "outs/agg/featurecounts_genes/countsstrandnonspecific_messy.txt",
        countsstrandnonspecific = "outs/agg/featurecounts_genes/countsstrandnonspecific.txt",
        metafeaturecounts = "outs/agg/featurecounts_genes/metafeature.counts.txt"
    params: 
        gtf = config['annotation_genes'],
        featureCountsstrandparam = getFeatureCountsStrandParam()
    conda:
        "omics"
    threads: 4
    shell: 
        """
featureCounts -s {params.featureCountsstrandparam} -T {threads} -t exon -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
featureCounts -s {params.featureCountsstrandparam} -T {threads} -B -O -a {params.gtf} -o {output.countsstrandnonspecificmessy} {input.sortedSTARbams}
cut -f1,7- {output.countsstrandnonspecificmessy} | awk 'NR > 1' > {output.countsstrandnonspecific}
featureCounts -T {threads} -B -O -a {params.gtf} -o {output.metafeaturecounts} {input.sortedSTARbams}
        """



rule featurecounts_genesandrtes:
#for non-telocal based counts for RTEs
    input:
        sortedSTARbams = expand("outs/{sample}/star_output/{sample}.sorted.primary.bam", sample = samples),
        libtype = "qc/library_type.txt"
    output:
        countsmessy = "outs/agg/featurecounts_genesandrtes/counts_messy.txt",
        counts = "outs/agg/featurecounts_genesandrtes/counts.txt",
    params: 
        gtf = "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations4/hs1.repMask.refseq.sorted.gtf",
        featureCountsstrandparam = getFeatureCountsStrandParam()
    conda:
        "omics"
    resources:
        cpus_per_task =10,
        runtime = 3000,
        mem_mb = 32000,
    shell: 
        """
featureCounts -s {params.featureCountsstrandparam} -M -T 8 --primary --ignoreDup --largestOverlap -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
        """