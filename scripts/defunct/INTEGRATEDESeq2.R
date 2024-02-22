log <- file(snakemake@log[[1]], open = "wt")
sink(log)


library("DESeq2")
library("readr")
library("dplyr")
library("pheatmap")
library("ggplot2")
library("tibble")
library("RColorBrewer")
library("cowplot")
library("PCAtools")



save.image()
### inputs

counttypes <- snakemake@params[["counttypes"]]
for (counttype in counttypes) {
    print(counttype)
    cts <- read.delim(snakemake@input[[counttype]])
    rownames(cts) <- cts$Geneid
    cts <- select(cts, -Geneid)
    cnames <- colnames(cts)

    coldata <- read.csv(snakemake@params[["peptable"]])

    contrasts <- snakemake@params[["contrasts"]]

    levels <- snakemake@params[["levels"]]

    outputdir <- snakemake@params[["outputdir"]]
    ###

    condition <- coldata$condition
    colnames(cts) <- coldata$sample_name
    dds <- DESeqDataSetFromMatrix(
        countData = cts,
        colData = coldata,
        design = ~ batch + condition
    )
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    # this sets prol as the reference level since its first in the vector
    dds$condition <- factor(dds$condition, levels = levels)

    ####
    dds <- DESeq(dds)
    save.image()
    ####
    resultsNames(dds) # lists the coefficients

    ####
    vst <- assay(vst(dds))
    vst_full <- vst(dds)
    sampleDists <- dist(t(vst))
    counttablesizenormed <- counts(dds, normalized = T)

    ## PCA plots
    pcaObj <- pca(vst, metadata = colData(dds), removeVar = 0.1)

    screep <- screeplot(pcaObj, title = "") +
        theme_cowplot() +
        theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

    pdf(paste(outputdir, counttype, "plots", "screeplot.pdf", sep = "/"), width = 5, height = 6)
    print(screep)
    dev.off()

    loadingsp <- plotloadings(pcaObj,
        components = getComponents(pcaObj, seq_len(3)),
        rangeRetain = 0.045, labSize = 4
    ) +
        theme(legend.position = "none") +
        theme_cowplot() +
        theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

    pdf(paste(outputdir, counttype, "plots", "loadings.pdf", sep = "/"), width = 6, height = 6)
    print(loadingsp)
    dev.off()

    pcaplot <- biplot(pcaObj,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
        colby = "condition", shape = "batch", legendPosition = "right",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    )
    pcap <- pcaplot +
        theme_cowplot() +
        theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

    pdf(paste(outputdir, counttype, "plots", "pcaplot.pdf", sep = "/"), width = 7, height = 6)
    print(pcap)
    dev.off()

    legend <- get_legend(
        # create some space to the left of the legend
        pcap + theme(legend.box.margin = margin(0, 0, 0, 1), legend.position = "right")
    )

    prow <- plot_grid(screep + theme(legend.position = "none", panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA)),
        loadingsp + theme(legend.position = "none", panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA)),
        pcap + theme(legend.position = "none", panel.border = element_rect(color = "black", linetype = 1, linewidth = 0.5, fill = NA)),
        nrow = 1,
        rel_widths = c(1, 1, 1), labels = "AUTO",
        align = "vh",
        axis = "bt"
    ) + theme(legend.position = "none")
    p <- plot_grid(prow, legend, nrow = 1, rel_widths = c(3, 0.4))

    pdf(paste(outputdir, counttype, "plots", "PCAgrid.pdf", sep = "/"), width = 16, height = 6)
    print(p)
    dev.off()

    # pdf("eigencorr.pdf", width=10, height=8)
    # eigencorplot(p, metavars = c('condition', 'sample_name'))
    # dev.off()


    pcaplot_statelipse <- biplot(pcaObj,
        colby = "condition", colkey = c("PRO" = "blue", "ESEN" = "yellow", "LSEN" = "red"),
        # ellipse config
        ellipse = TRUE,
        ellipseType = "t",
        ellipseLevel = 0.95,
        ellipseFill = TRUE,
        ellipseAlpha = 1 / 4,
        ellipseLineSize = 1.0,
        hline = 0, vline = c(-25, 0, 25),
        legendPosition = "top", legendLabSize = 16, legendIconSize = 8.0
    ) +
        theme_cowplot() +
        theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))
    pdf(paste(outputdir, counttype, "plots", "pcaplot_statelipse.pdf", sep = "/"), width = 8, height = 8)
    print(pcaplot_statelipse)
    dev.off()

    # PCA with batch effect removed

    vst_batchcorrected <- limma::removeBatchEffect(vst, colData(dds)$batch)

    pcaObj_batchcorrected <- pca(vst_batchcorrected, metadata = colData(dds), removeVar = 0.1)

    pcaplot <- biplot(pcaObj_batchcorrected,
        showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
        colby = "condition", shape = "batch", legendPosition = "right",
        labSize = 5, pointSize = 5, sizeLoadingsNames = 5
    )
    pcap <- pcaplot +
        theme_cowplot() +
        theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

    pdf(paste(outputdir, counttype, "plots", "pcaplot_batchcorrected.pdf", sep = "/"), width = 7, height = 6)
    print(pcap)
    dev.off()

    ############

    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vst_full$condition, vst_full$type, sep = "-")
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

    pdf(paste(outputdir, counttype, "plots", "heatmapplot.pdf", sep = "/"), width = 10, height = 8)
    pheatmap::pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        col = colors
    )
    dev.off()
    ####

    for (contrast in contrasts) {
        res <- results(dds, name = contrast)
        # or to shrink log fold changes association with condition:
        resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")
        # Save the plot as a PDF file
        pdf(paste(outputdir, counttype, contrast, "deplot.pdf", sep = "/"), width = 10, height = 8)
        plotMA(res, ylim = c(-7, 7), main = paste0(contrast))
        dev.off()
        pdf(paste(outputdir, counttype, contrast, "deLFCplot.pdf", sep = "/"), width = 10, height = 8)
        plotMA(resLFC, ylim = c(-10, 10), main = paste0(contrast))
        dev.off()


        resOrdered <- res[order(res$pvalue), ]
        resSig <- subset(resOrdered, padj < 0.1)
        counttablesizenormed <- counts(dds, normalized = T)
        rld <- rlog(dds, blind = FALSE)

        write.csv(as.data.frame(resOrdered), file = paste(outputdir, counttype, contrast, "results.csv", sep = "/"))
        write.csv(as.data.frame(resSig), file = paste(outputdir, counttype, contrast, "resultsSig.csv", sep = "/"))
        write.csv(as.data.frame(counttablesizenormed), file = paste(outputdir, counttype, contrast, "counttablesizenormed.csv", sep = "/"))
        write.csv(as.data.frame(assay(rld)), file = paste(outputdir, counttype, contrast, "rlogcounts.csv", sep = "/"))
    }
}

x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
