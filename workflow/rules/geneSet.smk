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
        r_annotation_families = config["r_annotation_families"],
        r_annotation_regions = config["r_annotation_regions"],
        r_annotation_intactness = config["r_annotation_intactness"],
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