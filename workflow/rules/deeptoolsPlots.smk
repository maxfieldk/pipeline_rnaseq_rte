
rule deeptools_AggcoverageMatrix:
    input:
        coverageF = expand("outs/{sample}/star_output/{sample}_forward_BPM_cov.bw", sample = samples),
        coverageR = expand("outs/{sample}/star_output/{sample}_reverse_BPM_cov.bw", sample = samples)
    params:
        l1hsorf2intact = config["l1hsorf2intact"]
    conda:
        "deeptools"
    log: "logs/agg/deeptools_AggcoverageMatrix.log"
    output:
        matrixF = "outs/agg/matrix_l1hsorf2intact_coverageF.tab.gz",
        matrixFsub = "outs/agg/matrix_l1hsorf2intact_coverage_subsetF.tab.gz",
        matrixR = "outs/agg/matrix_l1hsorf2intact_coverageR.tab.gz",
        matrixRsub = "outs/agg/matrix_l1hsorf2intact_coverage_subsetR.tab.gz",
        merged = "outs/agg/matrix_merged.tab.gz",
    shell:
        """
computeMatrix reference-point \
 -S {input.coverageF} \
 -R {params.l1hsorf2intact} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrixF} 2> {log}

computeMatrix reference-point \
 -S {input.coverageR} \
 -R {params.l1hsorf2intact} \
 --referencePoint TSS \
 --upstream 3000 \
 --binSize 100 \
 --downstream 9000 \
 -out {output.matrixR} 2> {log}

computeMatrixOperations info -m {output.matrixF} >> {log}
computeMatrixOperations info -m {output.matrixR} >> {log}
computeMatrixOperations filterStrand -m {output.matrixF} -o {output.matrixFsub} --strand +
computeMatrixOperations filterStrand -m {output.matrixR} -o {output.matrixRsub} --strand -
computeMatrixOperations rbind -m {output.matrixFsub} {output.matrixRsub} -o {output.merged}
        """

rule deeptools_AggcoverageMatrixScaled:
    input:
        coverageF = expand("outs/{sample}/star_output/{sample}_forward_BPM_cov.bw", sample = samples),
        coverageR = expand("outs/{sample}/star_output/{sample}_reverse_BPM_cov.bw", sample = samples)
    params:
        l1hsorf2intact = config["l1hsorf2intact"]
    conda:
        "deeptools"
    log: "logs/agg/deeptools_AggcoverageMatrixScaled.log"
    output:
        matrixF = "outs/agg/matrix_l1hsorf2intact_coverageFScaled.tab.gz",
        matrixFsub = "outs/agg/matrix_l1hsorf2intact_coverage_subsetFScaled.tab.gz",
        matrixR = "outs/agg/matrix_l1hsorf2intact_coverageRScaled.tab.gz",
        matrixRsub = "outs/agg/matrix_l1hsorf2intact_coverage_subsetRScaled.tab.gz",
        merged = "outs/agg/matrix_mergedScaled.tab.gz",
    shell:
        """
computeMatrix scale-regions \
 -S {input.coverageF} \
 -R {params.l1hsorf2intact} \
 -m 6000 \
 --upstream 0 \
 --downstream 0 \
 --binSize 10 \
 -out {output.matrixF} 2> {log}

computeMatrix scale-regions \
 -S {input.coverageR} \
 -R {params.l1hsorf2intact} \
 -m 6000 \
 --upstream 0 \
 --downstream 0 \
 --binSize 10 \
 -out {output.matrixR} 2> {log}

computeMatrixOperations info -m {output.matrixF} >> {log}
computeMatrixOperations info -m {output.matrixR} >> {log}
computeMatrixOperations filterStrand -m {output.matrixF} -o {output.matrixFsub} --strand +
computeMatrixOperations filterStrand -m {output.matrixR} -o {output.matrixRsub} --strand -
computeMatrixOperations rbind -m {output.matrixFsub} {output.matrixRsub} -o {output.merged}
        """

rule deeptools_plotAggcoverage:
    input:
        merged = "outs/agg/matrix_merged.tab.gz",
    conda:
        "deeptools"
    log: "logs/agg/deeptools_plotAggcoverage.log"
    output:
        merged = "results/agg/plots/heatmap_coverage_merged.png",

    shell:
        """
plotHeatmap \
 -m {input.merged} \
 -out {output.merged} \
 --heatmapHeight 15  \
 --refPointLabel TSS \
 --sortRegions descend \
 --sortUsing sum \
 --missingDataColor 0 \
 --labelRotation 0 \
 --linesAtTickMarks \
 --dpi 400 \
 --colorMap Blues \
 --samplesLabel PRO1 PRO2 PRO3 SEN1 SEN2 SEN3 \
 --yAxisLabel "Bins Per Million reads (BPM)" \
 --plotTitle "L1HS (ORF2 Intact)" 2>> {log}
        """

rule deeptools_plotAggcoverageScaled:
    input:
        merged = "outs/agg/matrix_mergedScaled.tab.gz",
    conda:
        "deeptools"
    log: "logs/agg/deeptools_plotAggcoverageScaled.log"
    output:
        merged = "results/agg/plots/heatmap_coverage_merged_scaled.png",

    shell:
        """
plotHeatmap \
 -m {input.merged} \
 -out {output.merged} \
 --heatmapHeight 15  \
 --sortRegions descend \
 --sortUsing sum \
 --missingDataColor 0 \
 --labelRotation 0 \
 --linesAtTickMarks \
 --dpi 400 \
 --colorMap Blues \
  --yAxisLabel "Bins Per Million reads (BPM)" \
 --samplesLabel PRO1 PRO2 PRO3 SEN1 SEN2 SEN3 \
 --plotTitle "L1HS (ORF2 Intact)" 2>> {log}
        """

rule deeptools_plotAggProfile:
    input:
        merged = "outs/agg/matrix_merged.tab.gz",
    conda:
        "deeptools"
    log: "logs/agg/deeptools_plotAggcoverage2.log"
    output:
        heatmap_coverage = "results/agg/plots/profile_coverage.png"
    shell:
        """
plotProfile -m {input.merged} \
--perGroup \
 --yAxisLabel "Bins Per Million reads (BPM)" \
--samplesLabel PRO1 PRO2 PRO3 SEN1 SEN2 SEN3 \
-out {output.heatmap_coverage} \
--dpi 400 \
--colors blue blue blue red red red \
--plotTitle "L1HS (ORF2 Intact)" 2>> {log}
        """

rule deeptools_plotAggProfileScaled:
    input:
        merged = "outs/agg/matrix_mergedScaled.tab.gz",
    conda:
        "deeptools"
    log: "logs/agg/deeptools_plotAggProfileScaled.log"
    output:
        heatmap_coverage = "results/agg/plots/profile_coverage_scaled.png"
    shell:
        """
plotProfile -m {input.merged} \
--perGroup \
 --yAxisLabel "Bins Per Million reads (BPM)" \
 --plotTitle "L1HS (ORF2 Intact)" \
-out {output.heatmap_coverage} \
--colors blue blue blue red red red \
--samplesLabel PRO1 PRO2 PRO3 SEN1 SEN2 SEN3 \
--dpi 400 2>> {log}
        """