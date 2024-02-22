source("~/data/common/myDefaults.r")

library("DESeq2")
library("readr")
library("pheatmap")
library("ggplot2")
library("tibble")
library("genefilter")
library("RColorBrewer")
library("cowplot")
library("PCAtools")
library("GGally")
library("tidyr")
library("bcbioRNASeq")
library("DESeqAnalysis")
library(AnnotationHub)
library(KEGGREST)
library(clusterProfiler)
library(AnnotationDbi)
library("biomaRt")
library(stringr)
library("dplyr")
library(EnhancedVolcano)

### functions
plotSave <- function(path, plot, width = 6, height = 6) {
    dir.create(dirname(path), recursive = TRUE)
    png(path, width = width, height = height, units = "in", res = 300)
    print(plot)
    dev.off()
}


conf <- c(
    confPrivate <- configr::read.config(file = "conf/private/configPrivate.yaml"),
    confShared <- configr::read.config(file = "conf/shared/configShared.yaml")
)


tryCatch(
    {
        params <- snakemake@params
        inputs <- snakemake@input
        outputs <- snakemake@output
    },
    error = function(e) {
        assign("params", list(
            "sample_table" = conf$sample_table,
            "tecounttype" = "telescope_multi",
            "contrasts" = conf$contrasts,
            "levels" = conf$levels,
            "outputdir" = "results/agg/deseq_telescope",
            "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
            "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation,
            "paralellize_bioc" = 8
        ), env = globalenv())

        assign("inputs", list(
            counts = "outs/agg/featurecounts_genes/counts.txt",
            rte_counts = sprintf("outs/%s/telescope/telescope-run_stats.tsv", read_delim(params$sample_table)$sample_name)
        ), env = globalenv())
    }
)

tecounttype <- params[["tecounttype"]]
print(tecounttype)
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]
coldata <- read.csv(params[["sample_table"]])
samples <- coldata$sample_name

if (params$paralellize_bioc) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

# Function for calculating size factors
# calculateSizeFactors <- function(mapped_frags) {
#     geomean <- expm1(mean(log1p(mapped_frags)))
#     mapped_frags / geomean
# }



# lib_size <- list()
# for (sample in samples) {
#     print(sample)
#     f <- sprintf("outs/%s/telescope/telescope-run_stats.tsv", sample)
#     h <- readLines(f, 1) %>%
#         strsplit(., "\t") %>%
#         unlist()
#     rstr <- sapply(strsplit(h[-c(1, 2)], ":"), function(t) as.numeric(unlist(t[2][1])))
#     names(rstr) <- sapply(strsplit(h[-c(1, 2)], ":"), function(t) t[1])
#     if (conf$READ_TYPE == "SE") {
#         rstr[["single_mapped"]]
#         lib_size[[sample]] <- rstr[["single_mapped"]]
#     } else {
#         rstr[["pair_mapped"]]
#         lib_size[[sample]] <- rstr[["pair_mapped"]]
#     }
# }

df <- read_delim(inputs[["rte_counts"]][1], comment = "#", col_names = FALSE)

if (tecounttype == "telescope_multi") {
    bounddf <- tibble(df[, 1]) %>% rename(gene_id = X1)
    for (sample in coldata$sample_name) {
        bounddf <- full_join(bounddf, read_delim(sprintf("outs/%s/telescope/telescope-run_stats.tsv", sample), comment = "#", col_names = FALSE) %>% select(X1, X3) %>% rename(gene_id = X1), by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", coldata$sample_name)
}

if (tecounttype == "telescope_unique") {
    bounddf <- tibble(df[, 1]) %>% rename(gene_id = X1)
    for (sample in coldata$sample_name) {
        bounddf <- full_join(bounddf, read_delim(sprintf("outs/%s/telescope/telescope-run_stats.tsv", sample), comment = "#", col_names = FALSE) %>% select(X1, X6) %>% rename(gene_id = X1), by = "gene_id")
    }
    colnames(bounddf) <- c("gene_id", coldata$sample_name)
}

bounddf1 <- bounddf[bounddf$gene_id != "__no_feature", ]

gene_cts <- read.delim(inputs$counts)
colnames(gene_cts) <- c("gene_id", coldata$sample_name)

cts <- rbind(gene_cts, as.data.frame(bounddf1 %>% replace(is.na(.), 0)))
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

# keep only genes with counts in at least one sample
# cts <- cts[rowSums(cts > 0) != 0, ]
# rounding since genes are allowed fractional counts
cts <- cts %>% mutate(across(everything(), ~ as.integer(round(.))))

condition <- coldata$condition
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~condition
)
# sizeFactors(dds) <- calculateSizeFactors(unlist(lib_size))

# I estimate the size factors using genes, and not RTEs, since there are 5M repeats and most have very very low counts
dds <- estimateSizeFactors(dds, controlGenes = rownames(dds) %in% gene_cts$gene_id)
# sizeFactors(dds)
# used to filter here
# this sets prol as the reference level since its first in the vector
dds$condition <- factor(dds$condition, levels = levels)

####
if (params$paralellize_bioc) {
    dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(params$paralellize_bioc))
} else {
    dds <- DESeq(dds)
}

keep <- rowSums(counts(dds)) >= 10
ddsfiltered <- dds[keep, ]
####
resultsNames(dds) # lists the coefficients

counttablesizenormed <- counts(dds, normalized = TRUE)

write.csv(as.data.frame(counttablesizenormed), file = paste(outputdir, tecounttype, "counttablesizenormed.csv", sep = "/"))
# write.csv(as.data.frame(assay(vst_assaydf)), file = paste(outputdir, tecounttype, "vstcounts.csv", sep = "/"))

# tag PLOTS
deseq_plots <- list()
for (contrast in contrasts) {
    res <- results(dds, name = contrast)
    res <- res[order(res$pvalue), ]
    write.csv(as.data.frame(res), file = paste(outputdir, tecounttype, contrast, "results.csv", sep = "/"))

    p <- EnhancedVolcano(res,
        lab = rownames(res),
        selectLab = c(""),
        title = contrast,
        drawConnectors = TRUE,
        x = "log2FoldChange",
        y = "padj",
        ylab = expression(-Log[10] ~ P["adj"]),
        legendPosition = "none",
        pCutoff = 0.05,
        xlim = c(-10, 10),
        ylim = c(-1, 15)
    ) + theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))
    mysave(paste(outputdir, tecounttype, contrast, "deplot.png", sep = "/"), 8, 8)
    deseq_plots[[tecounttype]][["volcano"]][[contrast]] <- p
}

vst <- vst(dds, blind = FALSE)
vst_assay <- assay(vst)

sampleDists <- dist(t(vst_assay))

## PCA plots
pcaObj <- pca(vst_assay, metadata = colData(dds), removeVar = 0.1)

p <- screeplot(pcaObj, title = "") +
    theme_cowplot() +
    mytheme
mysave(paste(outputdir, tecounttype, "screeplot.png", sep = "/"), 4, 4)
deseq_plots[[tecounttype]][["scree"]] <- p


p <- plotloadings(pcaObj,
    components = getComponents(pcaObj, seq_len(3)),
    rangeRetain = 0.045, labSize = 4
) +
    theme(legend.position = "none") +
    mytheme
mysave(paste(outputdir, tecounttype, "loadings.png", sep = "/"), 4, 4)
deseq_plots[[tecounttype]][["loadings"]] <- p


p <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
) +
    theme_gray() +
    mytheme
mysave(paste(outputdir, tecounttype, "pca.png", sep = "/"), 4, 4)
deseq_plots[[tecounttype]][["pca"]] <- p


p <- biplot(pcaObj,
    x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
) + mytheme
mysave(paste(outputdir, tecounttype, "pca34.png", sep = "/"), 4, 4)
deseq_plots[[tecounttype]][["pca34"]] <- p



sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$condition, vst$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

p <- pheatmap::pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
)
mysave(paste(outputdir, tecounttype, "pheatmap.png", sep = "/"), 4, 4)
deseq_plots[[tecounttype]][["dist_heatmap"]] <- p

save(deseq_plots, file = paste(outputdir, tecounttype, "deseq_plots.RData", sep = "/"))
save(dds, file = paste(outputdir, tecounttype, "dds.RData", sep = "/"))
