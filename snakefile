import os
import pandas as pd
from pathlib import Path

samples = config["samples"]
telocaltypes = config["telocaltypes"]
contrasts = config["contrasts"]
counttypes = config["counttypes"]
peptable = pd.read_csv("conf/private/peptable.csv")

#tag FILESTRUCTURE
paths = [
    "qc",
    "qc/{sample}",
    "rawdata/{sample}",
    "logs/{sample}",
    "outs/{sample}/trimmedReads",
    "outs/{sample}/{counttype}",
    "outs/agg/{counttype}",
    "results/{sample}/plots",
    "results/agg/clusterprofiler/{contrast}",
    "results/agg/deseq2/{counttype}/{contrast}",
    "results/agg/deseq_telescope/{tecounttype}/{contrast}",
    "results/agg/repeatanalysis_telescope"
    ]
paths = paths + []
for path in paths:
    for sample in samples:
        for counttype in config["counttypes"]:
            for tecounttype in config["tecounttypes"]:
                for contrast in config["contrasts"]:
                    os.makedirs(path.format(sample=sample, counttype=counttype,tecounttype=tecounttype, contrast=contrast), exist_ok = True)
#tag FUNCTIONS
def inferLibraryType():
        try: 
                with open("qc/library_type.txt") as f:
                        data = f.read()
                lines = re.findall('Fraction.*', data)
                pctMapped = [float(line.split(": ")[1]) for line in lines]
                if pctMapped[1] > 0.75:
                        libraryType = "forward"
                elif pctMapped[2] > 0.75:
                        libraryType = "reverse"
                else:
                        libraryType = "unstranded"
                return libraryType
        except:
                return "didn't run infer library type yet"

def getFeatureCountsStrandParam():
        libraryType = inferLibraryType()
        if libraryType == "forward":
                strandParam = "1"
        elif libraryType == "reverse":
                strandParam = "2"
        else:
                strandParam = "0"
        return strandParam

def getTElocalStrandParam():
        libraryType = inferLibraryType()
        if libraryType == "forward":
                strandParam = "forward"
        elif libraryType == "reverse":
                strandParam = "reverse"
        else:
                strandParam = "no"
        return strandParam

def getTelescopeStrandParam():
        libraryType = inferLibraryType()
        if libraryType == "forward":
                strandParam = "FR"
        elif libraryType == "reverse":
                strandParam = "RF"
        else:
                strandParam = "None"
        return strandParam


#tag PREPROCESSING
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

rule fastp_PE:
    input:
        r1=lambda wildcards: peptable.loc[peptable["sample_name"] == wildcards.sample, "file_path_R1"].iloc[0],
        r2=lambda wildcards: peptable.loc[peptable["sample_name"] == wildcards.sample, "file_path_R1"].iloc[0]
    priority: 100
    threads: 6
    log: "logs/{sample}/fastp.log"
    output:
        r1 = "outs/{sample}/trimmedReads/{sample}_1.trimmed.fastq.gz",
        r2 = "outs/{sample}/trimmedReads/{sample}_2.trimmed.fastq.gz",
        json = "outs/{sample}/trimmedReads/fastp.json",
        html = "outs/{sample}/trimmedReads/fastp.html"
    conda:
        "qc"
    shell:
        """
fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --json {output.json} --html {output.html} --detect_adapter_for_pe --thread {threads}
        """

rule fastp_SE:
    input:
        r1=lambda wildcards: peptable.loc[peptable["sample_name"] == wildcards.sample, "file_path_R1"].iloc[0],
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


#tag QC

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

#tag ALIGNMENT
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

#tag COUNTS
rule featurecounts_genes_PE:
    input:
        primaryBams = expand("outs/{sample}/star_output/{sample}.sorted.primary.bam", sample = samples),
        libtype = "qc/library_type.txt",
    output:
        countsmessy = "outs/agg/featurecounts_genes/counts_messy.txt",
        counts = "outs/agg/featurecounts_genes/counts.txt",
        readassignment = expand("outs/agg/featurecounts_genes/{sample}.sorted.primary.bam.featureCounts", sample = samples)

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

# rule filtertogetbamofnofeaturealignments:
#     input:
#         readassignment = "outs/agg/featurecounts_genes/{sample}.sorted.primary.bam.featureCounts",
#         sortedBam = "outs/{sample}/star_output/{sample}.sorted.bam"
#     output:
#         nofeaturefilteredbam = "outs/agg/featurecounts_genes/{sample}.sorted.nofeaturefiltered.bam"
#     conda:
#         "omics"
#     shell:
#         """
# awk '$2 ~ /Unassigned_NoFeatures/ {{print $1}}' {input.readassignment} > {input.readassignment}.nofeatureIDs.txt
# samtools view -b -N {input.readassignment}.nofeatureIDs.txt -o {output.nofeaturefilteredbam} {input.sortedBam}
# samtools stats {output.nofeaturefilteredbam} > {output.nofeaturefilteredbam}.stats
#         """


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

rule alignSTAR_SE:
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



rule featurecounts_genes_SE:
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



rule featurecounts_genesandrtes_SE:
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

#tag ALIGNMENT UTILITIES
rule common_sortIndexBam:
    input:
        sam = "outs/{sample}/star_output/Aligned.out.sam"
    output:
        sortedbam =  "outs/{sample}/star_output/{sample}.sorted.bam",
        stats =  "outs/{sample}/star_output/{sample}.sorted.bam.stats.txt",
        index = "outs/{sample}/star_output/{sample}.sorted.bam.bai"
    resources:
        cpus_per_task =10,
        mem_mb = 64000
    conda:
        "omics"
    shell:
        """
samtools view -@8 -b {input.sam} | \
samtools sort -@8 -m4g - > {output.sortedbam}
samtools index  -@6 {output.sortedbam}
samtools stats {output.sortedbam} > {output.stats}
        """


rule common_filterForPrimaryAlignments:
    input:
        "{path}.sorted.bam"
    output:
        bam =  "{path}.sorted.primary.bam",
        bamindex = "{path}.sorted.primary.bam.bai"
    threads: 4
    conda:
        "omics"
    shell: 
        """
samtools view -b -F 0x800 -F 0x100 -F 0x400 {input} > {output.bam}
samtools index {output.bam}
        """ 

#tag DIFFERENTIAL EXPRESSION

rule deseq:
    input:
        counts = "outs/agg/{counttype}/counts.txt",
    params:
        sample_table = config["sample_table"],
        contrasts = config["contrasts"],
        levels = config["levels"],
        paralellize_bioc = config["paralellize_bioc"],
        counttype = lambda w: w.counttype,
        outputdir = "results/agg/deseq2"
    resources:
        cpus_per_task =10,
        mem_mb = 200000,
        runtime = 1000
    conda: "deseq"
    output:
        results = expand("results/agg/deseq2/{{counttype}}/{contrast}/{resulttype}.csv", contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed"]),
        outfile = "results/agg/deseq2/{counttype}outfile.txt"
    script:
        "scripts/DESeq2.R"


rule deseq_telescope:
    input:
        counts = "outs/agg/featurecounts_genes/counts.txt",
        rte_counts = expand("outs/{sample}/telescope/telescope-run_stats.tsv", sample = samples)
    params:
        sample_table = config["sample_table"],
        contrasts = config["contrasts"],
        levels = config["levels"],
        paralellize_bioc = config["paralellize_bioc"],
        tecounttype = lambda w: w.tecounttype,
        outputdir = "results/agg/deseq_telescope"
    resources:
        cpus_per_task =10,
        mem_mb = 200000,
        runtime = 1000
    conda: "deseq"
    wildcard_constraints:
        tecounttype="[A-Za-z0-9_]+"
    output:
        results = expand("results/agg/deseq_telescope/{{tecounttype}}/{contrast}/results.csv", contrast = config["contrasts"]),
        counts_normed = "results/agg/deseq_telescope/{tecounttype}/counttablesizenormed.csv"
    script:
        "scripts/deseq_telescope.R"

#tag ENRICHMENT ANALYSIS
import os
rule enrichment_analysis:
    input:
        deresults = expand("results/agg/deseq2/featurecounts_genes/{contrast}/results.csv", contrast = config["contrasts"]),
        normcounttable = expand("results/agg/deseq2/featurecounts_genes/{contrast}/counttablesizenormed.csv", contrast = config["contrasts"])
    params:
        inputdir =lambda w, input: os.path.dirname(os.path.dirname(input.deresults[0])),
        contrasts = config["contrasts"],
        genesets_for_heatmaps = config["genesets_for_heatmaps"],
        genesets_for_gsea = config["genesets_for_gsea"],
        sample_table = config["sample_table"],
        outputdir = "results/agg/enrichment_analysis"
    conda:
        "ea"
    resources:
        cpus_per_task =10,
        mem_mb = 164000,
        runtime = 300
    output:
        outfile = "results/agg/enrichment_analysis/outfile.txt"
    script:
        "scripts/ea.R"

rule enrichment_analysis_repeats:
    input:
        resultsdf = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
    params:
        inputdir =lambda w, input: os.path.dirname(os.path.dirname(input.resultsdf[0])),
        r_annotation_fragmentsjoined = config["r_annotation_fragmentsjoined"],
        r_repeatmasker_annotation = config["r_repeatmasker_annotation"],
        contrasts = config["contrasts"],
        tecounttypes = config["tecounttypes"],
        sample_table = config["sample_table"],
        outputdir = lambda w, output: os.path.dirname(output.outfile)
    conda:
        "ea"
    resources:
        cpus_per_task =10,
        mem_mb = 164000,
        runtime = 300
    output:
        outfile = "results/agg/enrichment_analysis_repeats/{tecounttype}/outfile.txt"
    script:
        "scripts/ea_repeats.R"

#tag REPETITIVE ELEMENTS
rule collateBam:
    input:
        sortedBam = "outs/{sample}/star_output/{sample}.sorted.bam"
    output:
        collatedbam = "outs/{sample}/star_output/{sample}.collated.bam",
    conda:
        "omics"
    resources:
        mem_mb  = 128000,
        runtime = 60
    shell:
        """
samtools collate -o {output.collatedbam} {input.sortedBam}
        """

rule telescope:
    input:
        collatedbam = "outs/{sample}/star_output/{sample}.collated.bam",
        libtype = "qc/library_type.txt"
    params:
        gtf = config["annotation_rtes"],
        strandparam = getTelescopeStrandParam(),
    output:
        counts = "outs/{sample}/telescope/telescope-run_stats.tsv"
    threads: 4
    conda:
        "telescope3"
    resources:
        mem_mb  = 128000,
        runtime = 600
    shell: 
        """
telescope assign \
--attribute gene_id \
--ncpu 1 \
--stranded_mode {params.strandparam} \
--outdir $(dirname {output.counts}) \
{input.collatedbam} \
{params.gtf}
        """

rule repeatanalysis_telescope:
    input:
        results = expand("results/agg/deseq_telescope/{tecounttype}/{contrast}/results.csv", contrast = config["contrasts"], tecounttype = config["tecounttypes"]),
        counts_normed = expand("results/agg/deseq_telescope/{tecounttype}/counttablesizenormed.csv", tecounttype = config["tecounttypes"])
    params:
        inputdir = "results/agg/deseq_telescope",
        outputdir = "results/agg/repeatanalysis_telescope"
    conda:
        "repeatanalysis"
    resources:
        cpus_per_task =10,
        mem_mb = 164000,
        runtime = 300
    output:
        resultsdf = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
    script:
        "scripts/repeatanalysis.R"

rule repeatanalysis_plots:
    input:
        resultsdf = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
    params:
        r_annotation_fragmentsjoined = config["r_annotation_fragmentsjoined"],
        repeatmasker_annotation = config["r_repeatmasker_annotation"],
        contrasts = config["contrasts"],
        tecounttypes = config["tecounttypes"],
        levelslegendmap = config["levelslegendmap"],
        inputdir = "results/agg/deseq_telescope",
        outputdir = "results/agg/repeatanalysis_telescope"
    conda:
        "repeatanalysis"
    resources:
        cpus_per_task =10,
        mem_mb = 164000,
        runtime = 300
    output:
        outfile = "results/agg/repeatanalysis/plots.outfile.txt"
    script:
        "scripts/repeatanalysisPlots.R"
    

rule repeatVariance:
    input: 
        resultsdf = "results/agg/repeatanalysis/resultsdf.tsv"
    params:
        repeatanalysis = config["repeatanalysis"],
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        telocaltypes = config["telocaltypes"],
        levelslegendmap = config["levelslegendmap"],
        lengthreq = config["lengthreq"],
        peptable = "conf/private/peptable.csv"
    output:
        corrplot1 = "results/agg/repeatanalysis/variance/{telocaltype}/{contrast}/corrplot1.pdf",
        corrplot2 = "results/agg/repeatanalysis/variance/{telocaltype}/{contrast}/corrplot2.pdf",
        corrplotcontrast = "results/agg/repeatanalysis/variance/{telocaltype}/{contrast}/corrplot{contrast}.pdf"
    notebook:
        "scripts/rteVariance.ipynb"