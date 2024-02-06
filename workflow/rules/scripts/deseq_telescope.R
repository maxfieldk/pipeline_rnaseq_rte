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
            "paralellize_bioc" = 8
        ), env = globalenv())

        assign("inputs", list(
            counts = sprintf("outs/%s/telescope/telescope-run_stats.tsv", read_delim(params$sample_table)$sample_name)
        ), env = globalenv())
    }
)

tecounttype <- params[["tecounttype"]]
contrasts <- params[["contrasts"]]
levels <- params[["levels"]]
outputdir <- params[["outputdir"]]


if (params$paralellize_bioc) {
    library(BiocParallel)
    register(MulticoreParam(8))
}

# Function for calculating size factors
calculateSizeFactors <- function(mapped_frags) {
    geomean <- expm1(mean(log1p(mapped_frags)))
    mapped_frags / geomean
}

coldata <- read.csv(params[["sample_table"]])
samples <- coldata$sample_name

lib_size <- list()
for (sample in samples) {
    print(sample)
    f <- sprintf("outs/%s/telescope/telescope-run_stats.tsv", sample)
    h <- readLines(f, 1) %>%
        strsplit(., "\t") %>%
        unlist()
    rstr <- sapply(strsplit(h[-c(1, 2)], ":"), function(t) as.numeric(unlist(t[2][1])))
    names(rstr) <- sapply(strsplit(h[-c(1, 2)], ":"), function(t) t[1])
    if (conf$READ_TYPE == "SE") {
        rstr[["single_mapped"]]
        lib_size[[sample]] <- rstr[["single_mapped"]]
    } else {
        rstr[["pair_mapped"]]
        lib_size[[sample]] <- rstr[["pair_mapped"]]
    }
}

df <- read_delim(inputs[["counts"]][1], comment = "#", col_names = FALSE)
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


cts <- as.data.frame(bounddf %>% replace(is.na(.), 0))
rownames(cts) <- cts$gene_id
cts <- dplyr::select(cts, -gene_id)
cnames <- colnames(cts)

condition <- coldata$condition
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~condition
)
sizeFactors(dds) <- calculateSizeFactors(unlist(lib_size))

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
counttablesizenormed <- counts(dds, normalized = T)

vst_assay <- vst(dds, blind = FALSE)
write.csv(as.data.frame(counttablesizenormed), file = paste(outputdir, tecounttype, "counttablesizenormed.csv", sep = "/"))
write.csv(as.data.frame(assay(vst_assay)), file = paste(outputdir, tecounttype, "vstcounts.csv", sep = "/"))

############
for (contrast in contrasts) {
    res <- results(dds, name = contrast)
    res <- res[order(res$pvalue), ]
    write.csv(as.data.frame(res), file = paste(outputdir, tecounttype, contrast, "results.csv", sep = "/"))
}
