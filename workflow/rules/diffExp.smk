
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
        counts = expand("outs/{sample}/telescope/telescope-run_stats.tsv", sample = samples)
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
    output:
        results = expand("results/agg/deseq_telescope/{{tecounttype}}/{contrast}/{resulttype}.csv", contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed"])
    script:
        "scripts/deseq_telescope.R"