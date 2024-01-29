
rule collateBam:
    input:
        sortedSTARbam = "outs/{sample}/star_output/{sample}.sorted.bam",
    output:
        collatedbam = "outs/{sample}/star_output/{sample}.collated.bam",
    conda:
        "omics"
    shell:
        """
samtools collate -o {output.collatedbam} {input.sortedSTARbam}
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

# rule TElocal:
#     input:
#         sortedSTARbam = "outs/{sample}/star_output/{sample}.sorted.bam",
#         sortedbamSTARindex = "outs/{sample}/star_output/{sample}.sorted.bam.bai",
#         libtype = "qc/library_type.txt"

#     log: "logs/{sample}/TElocal_{maptype}.log"
#     params:
#         annotation_genes = config["annotation_genes"],
#         locindTElocal = config["locindTElocal"],
#         telocalstrandparam = getTElocalStrandParam(),
#         outputprefix = "outs/{sample}/telocal_{maptype}/{sample}",
#     output:
#         counts = "outs/{sample}/telocal_{maptype}/{sample}.cntTable",
#     threads: 4
#     conda:
#         "TElocal"
#     resources:
#         mem_mb  = 40000
#     shell: 
#         """
# TElocal --sortByPos -b {input.sortedSTARbam} --stranded {params.telocalstrandparam} --mode {wildcards.maptype} --GTF {params.annotation_genes} --TE {params.locindTElocal} --project {params.outputprefix}
#         """

# rule mergeTElocal:
#     input:
#         telocal = expand("outs/{sample}/telocal_{{maptype}}/{sample}.cntTable", sample = samples)
#     params:
#         sample_table = config["sample_table"]
#     conda: "renv"
#     output:
#         aggcounts = "outs/agg/telocal_{maptype}/counts.txt",
#     script:
#         "scripts/mergeTEdf.R"



# rule mergeTElocal:
#     input:
#         telocal = expand("outs/{sample}/TElocal/{sample}_{{maptype}}.cntTable", sample = samples)
#     params:
#         sample_table = config["sample_table"]
#     conda: "renv"
#     output:
#         aggcounts = "outs/agg/TElocalCounts_{maptype}.txt",
#     shell:
#         """
# samples=({input.telocal})
# cut -f 1 ${{samples[1]}} > combined_output.tsv
# for file in ${{samples[@]}}
# do
#   # Extract the second column and append to the growing output file
#   cut -f 2 "$file" | paste combined_output.tsv - > temp_output.tsv
#   mv temp_output.tsv combined_output.tsv
# done
# mv combined_output.tsv {output.aggcounts}
#         """


rule repeatanalysis:
    input:
        deseq = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
    params:
        inputdir = "results/agg/deseq",
        outputdir = "results/agg/repeatanalysis"
    conda:
        "repeatanalysis"
    resources:
        cpus_per_task =10,
        mem_mb = 164000,
        runtime = 300
    log:
        "logs/agg/repeatanalysis.log"
    output:
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv",
        resultsdf = "results/agg/repeatanalysis/resultsdf.tsv",
        sigAluYs = "results/agg/repeatanalysis/sigAluYs.tsv",
        sigL1s = "results/agg/repeatanalysis/sigL1s.tsv",
        sigHERVKs = "results/agg/repeatanalysis/sigHERVKs.tsv",
        outfile = "results/agg/repeatanalysis/outfile.txt"
    script:
        "scripts/repeatanalysis.R"


rule repeatanalysis_telescope:
    input:
        deseq = expand("results/agg/deseq_telescope/{tecounttypes}/{contrast}/results.csv", tecounttypes = config["tecounttypes"], contrast = config["contrasts"]),
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
        "evo2"
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