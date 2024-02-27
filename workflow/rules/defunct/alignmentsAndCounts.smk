rule alignSTAR_PE:
    input:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz"
    params:
        index = config["starindex"],
        outdir = "outs/{sample}/star_output/"
    output:
        bam = temp("outs/{sample}/star_output/Aligned.out.sam")
    log:
        out = "logs/{sample}/STAR.out",
        err = "logs/{sample}/STAR.err"
    threads: 8
    resources:
        mem_mb  = 60000
    conda:
        "star"
    shell:
        """
STAR --genomeDir {params.index} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {params.outdir} --runThreadN {threads} --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 > {log.out} 2> {log.err}
        """

rule featurecounts_genes_PE:
    input:
        primaryBams = expand("outs/{sample}/star_output/{sample}.sorted.primary.bam", sample = samples),
        libtype = "qc/library_type.txt",
    output:
        countsmessy = "outs/agg/featurecounts_genes/counts_messy.txt",
        counts = "outs/agg/featurecounts_genes/counts.txt",
        readassignment = expand("outs/{sample}/star_output/{sample}.sorted.primary.bam.featureCounts", sample = samples)

    params: 
        gtf = config['annotation_genes'],
        featureCountsstrandparam = getFeatureCountsStrandParam()
    conda: "omics"
    threads: 4
    shell: 
        """
featureCounts -p -s {params.featureCountsstrandparam} -O -T {threads} -t exon -a {params.gtf} -o {output.countsmessy} -R CORE  --minOverlap 25 --fraction --primary {input.primaryBams}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
        """

rule filtertogetbamofnofeaturealignments:
    input:
        readassignment = "outs/{sample}/star_output/{sample}.sorted.primary.bam.featureCounts",
        sortedBam = "outs/agg/featurecounts_genes/{sample}.sorted.bam"
    output:
        nofeaturefilteredbam = "outs/agg/featurecounts_genes/{sample}.sorted.nofeaturefiltered.bam"
    conda:
        "omics"
    shell:
        """
awk '$2 ~ /Unassigned_NoFeatures/ {{print $1}}' {input.readassignment} > {input.readassignment}.nofeatureIDs.txt
samtools view -b -N {input.readassignment}.nofeatureIDs.txt -o {output.nofeaturefilteredbam} {input.sortedBam}
samtools stats {output.nofeaturefilteredbam} > {output.nofeaturefilteredbam}.stats
        """


rule featurecounts_genesandrtes_PE:
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
featureCounts -p -s {params.featureCountsstrandparam} -M -T 8 --primary --ignoreDup --largestOverlap -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
        """