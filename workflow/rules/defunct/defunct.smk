#######################################################################
#Bowtie based alignment stats
#######################################################################
rule repeatCounts:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam",
        bamindex = "outs/{sample}/{sample}.sorted.bam.bai",
        repeats = config["repeats"]
    log: "logs/{sample}/repeatCounts.log"
    output:
        repeatcounts = "outs/{sample}/{sample}.repeatmasker.counts"
    threads: 4
    conda:
        "deeptools"
    shell:
        "bedtools coverage -counts -F 0.5 -sorted -b {input.bam} -a {input.repeats} > {output.repeatcounts} 2> {log}"


rule countMappednonmitoreads:
    input:
        bam = "outs/{sample}/{sample}.sorted.bam"
    output:
        mnma = "outs/{sample}/{sample}_number_nonmito_mapped_per_million.txt"
    log: "logs/{sample}/countMappednonmitoreads.log"
    conda:
        "omics"
    shell:
        "samtools view -F 4 {input.bam} | awk '$3 !~ /^MT/ { count++ } END { print count/1000000 }' > {output.mnma} "

rule getnormalizedCounts:
    input:
        repeatcounts = "outs/{sample}/{sample}.repeatmasker.counts",
        mnma = "outs/{sample}/{sample}_number_nonmito_mapped_per_million.txt"
    threads: 4
    output:
        repeatcounts = "results/{sample}/{sample}.repeatmasker.countsPerMillionNonMito"
    log: "logs/{sample}/getnormalizedCounts.log"
    conda:
        "deeptools"
    shell:
        """
awk -v milnonmito=$(cat {input.mnma}) -F "\t" '{{print $0 "\t" $NF/milnonmito}}' {input.repeatcounts} > {output.repeatcounts} 2> {log}
        """

rule getnormalizedFamilyCounts:
    input:
        counts = "results/{sample}/{sample}.repeatmasker.countsPerMillionNonMito"
    threads: 4
    output:
        l1 = "results/{sample}/{sample}.l1.counts.txt",
        l1hs = "results/{sample}/{sample}.l1hs.counts.txt",
        alu = "results/{sample}/{sample}.alu.counts.txt",
        aluy = "results/{sample}/{sample}.aluy.counts.txt",
        herv = "results/{sample}/{sample}.herv.counts.txt",
        hervk = "results/{sample}/{sample}.hervk.counts.txt",
        sva = "results/{sample}/{sample}.sva.counts.txt"
    conda:
        "deeptools"
    shell:
        """
        awk -F "\t" '/LINE\/L1/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.l1}
awk -F "\t" '/L1HS/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.l1hs}
awk -F "\t" '/SINE\/Alu/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.alu}
awk -F "\t" '/AluY/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.aluy}
awk -F "\t" '/LTR\/ERV/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.herv}
awk -F "\t" '/HERVK/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.hervk}
awk -F "\t" '/SVA/ {{sum+=$NF}} END {{print sum}}' {input.counts} > {output.sva}
        """

rule multiBamSummary:
    input:
        bam = expand("outs/{a}/{a}.sorted.bam", a=samples)
    output:
        agg = "results/agg/multiBamSummary.npz"
    threads: 4
    log: "logs/multiBamSummary.log"
    conda:
        "deeptools"
    shell: "multiBamSummary bins --numberOfProcessors {threads} --bamfiles {input.bam} -o {output.agg} 2> {log}"

rule plotbamPCA:
    input:
        agg = "results/agg/multiBamSummary.npz"
    output:
        plot = report("results/agg/plots/PCA_readCounts.png", caption = "report/plotbamPCA.rst", category="PCA")
    log: "logs/agg/plotbamPCA.log"
    conda:
        "deeptools"
    shell:
        "plotPCA -in {input.agg} -o {output.plot} -T 'PCA of read counts' 2> {log}"



rule gvizreduxParallel:
    input:
        sortedSTARbams = expand("outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = samples),
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv"
    params:
        contrasts = config["contrasts"],
        telocaltypes = config["telocaltypes"],
        peptable = "conf/private/peptable.csv",
        refseq = config["refseq2"],
        genes = config["genes"],
        ideogram = config["ideogram"],
        intactl1s = config["intactl1s"],
        repeats = config["repeats2"],
        levels = config["levels"],
        condition_colors = config["condition_colors"],
        telocalmapping = "/users/mkelsey/data/ref/genomes/hs1/TElocal/T2T_CHM13_v2_rmsk_TE.gtf.locInd.locations",
        outputdir = "results/agg/genometracks"
    resources:
        mem_mb  = 20000,
        runtime = 30
    output:
        outfile = "temp/outfile{telocaltype}{contrast}{rtekind}{range}._{telocaltype}._{contrast}._{rtekind}._{range}.gviz.txt"
    conda:
        "repeatanalysis"
    script:
        "scripts/gviz.R"



rule ideogram:
    input:
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv"
    params:
        contrasts = config["contrasts"],
        telocaltypes = config["telocaltypes"],
        rtestoplot = config["rtestoplot"],
        sample_table = config["sample_table"],
        karyotype = config["karyotype"],
        genedensity = config["genedensity"],
        namedcolorlist = config["namedcolorlist"],
        namedmarkerlist = config["namedmarkerlist"],
        outputdir = "results/agg/repeatanalysis"
    conda:
        "repeatanalysis"
    log:
        "logs/agg/ideogram.log"
    output:
        outfile = "results/agg/repeatanalysis/ideogram_outfile.txt"
    script:
        "scripts/ideogram.R"


rule getDEelementMSA:
    input:
        "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.bed"
    params:
        consensus = "/users/mkelsey/data/ref/sequences/{rte}consensus.fa",
        hs1fa = config["reference"],
        outgroup = "/users/mkelsey/data/ref/sequences/{rte}outgroup.fa"
    conda:
        "evo"
    output:
        fa = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.fasta",
        faWithConsensus = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.withconsensus.fasta",
        faWithConsensusAndOutgroup = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.withconsensusandoutgroup.fasta",
        aln = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.aln",
        alnWithConsensus = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.WithConsensus.aln",
        alnWithConsensusAndOutgroup = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{rte}/de{direction}.WithConsensusAndOutgroup.aln"

    shell:
        """
if [ -s {input} ]
then
bedtools getfasta -s -fullHeader -fi {params.hs1fa} -bed {input} -fo {output.fa}
cat {params.consensus} > {output.faWithConsensus}
echo >> {output.faWithConsensus}
cat {output.fa} >> {output.faWithConsensus}
cat {params.consensus} > {output.faWithConsensusAndOutgroup} 
echo >> {output.faWithConsensusAndOutgroup}
cat {params.outgroup} >> {output.faWithConsensusAndOutgroup}
echo >> {output.faWithConsensusAndOutgroup}
cat {output.fa} >> {output.faWithConsensusAndOutgroup}
if [ $(grep ">" {output.fa} | wc -l) -gt 1 ]
then
mafft --auto {output.fa} > {output.aln}
else
touch {output.aln}
fi
mafft --auto {output.faWithConsensus} > {output.alnWithConsensus}
mafft --auto {output.faWithConsensusAndOutgroup} > {output.alnWithConsensusAndOutgroup}
else
touch {output.fa}
touch {output.faWithConsensus}
touch {output.faWithConsensusAndOutgroup}
touch {output.aln}
touch {output.alnWithConsensus}
touch {output.alnWithConsensusAndOutgroup}
fi
        """ 
#expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.fasta", telocaltype = telocaltypes, contrast = contrasts, RTE = "L1HS", direciton = ["UP", "DOWN"])


rule evoAnalysis:
    input:
        alnWithConsensus = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.WithConsensus.aln",
        aln = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.aln"
    conda:
        "evo"
    params:
        basepath = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}"
    output:
        alignment = "results/agg/repeatanalysis/{telocaltype}/{contrast}/{RTE}/de{direction}.alignment.pdf"
    script:
        "scripts/evoAnalysis.R"




#####################################################################################For the integration project only
try:
    arnasamples = peptable[peptable.batch == "alexandra"].sample_name
    marco2samples = peptable[peptable.batch == "marco2"].sample_name
    nat2019samples = peptable[peptable.batch == "nat2019"].sample_name
except:
    arnasamples = "a"
    marco2samples = "a"
    nat2019samples = "a"

### BE SURE to order samples in the input according to their order in the sample_table.csv
rule INTEGRATEfeatureCounts:
    input:
        sortedSTARbams = expand("/users/mkelsey/data/arna/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = arnasamples) + expand("/users/mkelsey/data/senescence/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = marco2samples) + expand("/users/mkelsey/data/marco/outs/{sample}/star_output/{sample}.STAR.sorted.bam", sample = nat2019samples)
    output:
        countsmessy = "outs/agg/INTEGRATE_refseq.counts_messy.txt",
        counts = "outs/agg/INTEGRATE_refseq.counts.txt",
        metafeaturecounts = "outs/agg/INTEGRATE_refseq.metafeature.counts.txt"
    params: 
        gtf = config['refseq']
    log: "logs/agg/featureCounts.log"
    conda:
        "deeptools"
    threads: 2
    shell: 
        """
featureCounts -p -T {threads} -t exon -a {params.gtf} -o {output.countsmessy} {input.sortedSTARbams} 2> {log}
cut -f1,7- {output.countsmessy} | awk 'NR > 1' > {output.counts}
featureCounts -p -T {threads} -B -O -a {params.gtf} -o {output.metafeaturecounts} {input.sortedSTARbams} 2>> {log}
        """

### BE SURE to order samples in the input according to their order in the sample_table.csv
rule INTEGRATEmergeTElocal:
    input:
        telocal_multi = expand("/users/mkelsey/data/arna/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = arnasamples) + expand("/users/mkelsey/data/senescence/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = marco2samples) + expand("/users/mkelsey/data/marco/outs/{sample}/TElocal/{sample}_multi.cntTable", sample = nat2019samples),
        telocal_uniq =  expand("/users/mkelsey/data/arna/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = arnasamples) + expand("/users/mkelsey/data/senescence/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = marco2samples) + expand("/users/mkelsey/data/marco/outs/{sample}/TElocal/{sample}_uniq.cntTable", sample = nat2019samples)
    params:
        sample_table = config["sample_table"],
        telocaltypes = config["telocaltypes"]
    conda: "renv"
    log: "logs/agg/INTEGRATEmergeTElocal.log"
    output:
        telocal_multi = "outs/agg/INTEGRATE_TElocalCounts_multi.txt",
        telocal_uniq = "outs/agg/INTEGRATE_TElocalCounts_uniq.txt"
    script:
        "scripts/INTEGRATEmergeTEdf.R"    


rule INTEGRATEDEseq2:
    input:
        star = "outs/agg/INTEGRATE_refseq.counts.txt",
        telocal_multi = "outs/agg/INTEGRATE_TElocalCounts_multi.txt",
        telocal_uniq = "outs/agg/INTEGRATE_TElocalCounts_uniq.txt",
    params:
        peptable = "conf/private/peptable.csv",
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        levels = config["levels"],
        outputdir = "results/agg/deseq2"
    conda: "renv"
    log: "logs/agg/DEseq2.log"
    output:
        results = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
        outfile = "results/agg/deseq2/INTEGRATE_outfile.txt"
    script:
        "scripts/INTEGRATEDESeq2.R"


#be sure to order the contrasts if there are multiple by the one whose naming scheme is respected in this project
rule INTEGRATEdetermineSharedDeTes:
    params:
        peptable = "conf/private/peptable.csv",
        outputdir = "results/agg/deseq2",
        telocaltypes = config["telocaltypes"],
        contraststocompare = ["condition_SEN_vs_PRO", "condition_LSEN_vs_PRO"]
    conda: "repeatanalysis"
    output:
        outfile = "outfile_sharedetes.txt"
    script:
        "scripts/INTEGRATEanalyzeDERTEs.R"


rule INTEGRATErepeatanalysis:
    input:
        deseq = expand("results/agg/deseq2/{counttype}/{contrast}/{resulttype}.csv", counttype = config["counttypes"], contrast = config["contrasts"], resulttype = ["results", "counttablesizenormed", "rlogcounts"]),
        telocal = expand("outs/agg/INTEGRATE_TElocalCounts_{maptype}.txt", sample = samples, maptype = ["multi", "uniq"])
    params:
        contrasts = config["contrasts"],
        counttypes = config["counttypes"],
        telocaltypes = config["telocaltypes"],
        levelslegendmap = config["levelslegendmap"],
        peptable = "conf/private/peptable.csv",
        contrast_colors =config["contrast_colors"],
        condition_colors =config["condition_colors"],
        repeats = config["repeats"],
        telocalmapping = config["telocalmapping"],
        inputdir = "results/agg/deseq2",
        outputdir = "results/agg/repeatanalysis"
    conda:
        "repeatanalysis"
    log:
        "logs/agg/repeatanalysis.log"
    output:
        activeelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/activeelementContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisactiveelementContrastplot.rst", category="repeat analysis"),
        familyContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/FamilyContrastplot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysisfamilyContrastplot.rst", category="repeat analysis"),
        combinedelementContrastplot = report(expand("results/agg/repeatanalysis/{telocaltype}/{contrast}/CombinedContrastPlot.pdf", telocaltype = config["telocaltypes"], contrast = config["contrasts"]),caption = "report/repeatanalysiscombinedContrastplot.rst", category="repeat analysis"),
        sharedamongallcontrasts_derte = "results/agg/repeatanalysis/sharedamongallcontrasts_derte.tsv",
        DETEsbyContrast = "results/agg/repeatanalysis/allactiveDETEs.tsv",
        outfile = "results/agg/repeatanalysis/INTEGRATEoutfile.txt"
    script:
        "scripts/repeatanalysis.R"



#############################################

