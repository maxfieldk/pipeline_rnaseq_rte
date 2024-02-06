library("ggplot2")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("org.Hs.eg.db")
library("ggrepel")
library("grid")
library("readr")
library("stringr")
library("dplyr")
library("tibble")
library("tidyr")

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
            "inputdir" = "results/agg/deseq_telescope",
            "outputdir" = "results/agg/repeatanalysis_telescope"
        ), env = globalenv())
        assign("outputs", list(
            "resultsdf" = "results/agg/repeatanalysis_telescope/resultsdf.tsv",
            "vstresultsdf" = "results/agg/repeatanalysis_telescope/resultsdfvst.tsv"
        ), env = globalenv())
    }
)



outputdir <- params$outputdir
dir.create(outputdir, recursive = TRUE, showWarnings = FALSE)
# order matters for the colors!
contrast_colors <- conf$contrast_colors
condition_colors <- conf$condition_colors
contrasts <- conf$contrasts
levelslegendmap <- conf$levelslegendmap
tecounttypes <- conf$tecounttypes
tecounttypes <- conf$tecounttypes
lengthreq <- conf$lengthreq
maincontrast <- contrasts[1]



contrast_colors <- unname(unlist(contrast_colors))
condition_colors <- unname(unlist(condition_colors))
peptable <- read.csv(conf$peptable)




##########################



# build dds frames
RESLIST <- list()
for (tecounttype in tecounttypes) {
    contrastl <- list()
    for (contrast in contrasts) {
        ddsres <- read_csv(paste(params$inputdir, tecounttype, contrast, "results.csv", sep = "/"))
        ddsresmod <- ddsres %>%
            dplyr::rename(gene_id = ...1) %>%
            mutate(Significance = ifelse(padj < 0.05, ifelse(padj < 0.001, "< 0.001", "< 0.05"), "> 0.05")) %>%
            dplyr::select(c(gene_id, log2FoldChange, stat, padj, Significance)) %>%
            dplyr::rename(!!paste0("log2FoldChange_", contrast) := log2FoldChange) %>%
            dplyr::rename(!!paste0("stat_", contrast) := stat) %>%
            dplyr::rename(!!paste0("padj_", contrast) := padj) %>%
            dplyr::rename(!!paste0("Significance_", contrast) := Significance)
        contrastl[[contrast]] <- ddsresmod
    }
    ddsrestetype <- Reduce(function(x, y) merge(x, y, by = "gene_id"), contrastl, accumulate = FALSE)
    ddsrestetype <- ddsrestetype %>% mutate(tecounttype = tecounttype)
    RESLIST[[tecounttype]] <- ddsrestetype
}
ddsfinal <- bind_rows(RESLIST)

COUNTLIST <- list()
VSTCOUNTLIST <- list()
for (tecounttype in tecounttypes) {
    conditions <- peptable$condition %>% unique()
    ddscounts <- read_csv(paste(params$inputdir, tecounttype, "counttablesizenormed.csv", sep = "/"))
    ddsvstcounts <- read_csv(paste(params$inputdir, tecounttype, "vstcounts.csv", sep = "/"))
    colnames(ddscounts)[1] <- "gene_id"
    colnames(ddsvstcounts)[1] <- "gene_id"
    avoidzero <- 1
    meancols <- c()
    log2meancols <- c()
    vstmeancols <- c()
    for (condition in conditions) {
        condition_samples <- filter(peptable, condition == {{ condition }})$sample_name
        s <- ""
        for (e in condition_samples) {
            if (s == "") {
                s <- paste0("`", e, "`")
            } else {
                s <- paste0(s, "+", "`", e, "`")
            }
        }
        s <- paste0("(", s, ")", "/", length(condition_samples))
        meancol <- paste0(condition, "mean")
        ddscounts <- ddscounts %>% mutate({{ meancol }} := eval(parse(text = s)))
        vstmeancol <- paste0("vst", condition, "mean")
        ddsvstcounts <- ddsvstcounts %>% mutate({{ vstmeancol }} := eval(parse(text = s)))
        log2meancol <- paste0("log2", condition, "mean")
        ddscounts <- ddscounts %>% mutate({{ log2meancol }} := log2(.data[[meancol]] + avoidzero))
        meancols <- c(meancols, meancol)
        log2meancols <- c(log2meancols, log2meancol)
        vstmeancols <- c(vstmeancols, vstmeancol)
    }
    vstcolumns_to_retain <- c("gene_id", vstmeancols)
    log2sub <- ddscounts
    vstsub <- ddsvstcounts[, vstcolumns_to_retain]
    log2sub <- log2sub %>% mutate(tecounttype = tecounttype)
    countres <- merge(log2sub, vstsub, by = "gene_id")
    COUNTLIST[[tecounttype]] <- countres
    ddsvstcounts <- ddsvstcounts %>% mutate(tecounttype = tecounttype)
    VSTCOUNTLIST[[tecounttype]] <- ddsvstcounts
}
countsfinal <- bind_rows(COUNTLIST)
vstcountsfinal <- bind_rows(VSTCOUNTLIST)
# merge counts and dds, then add TE annotations


# merge counts and dds, then add TE annotations
resultsdf <- merge(ddsfinal, countsfinal, by = c("gene_id", "tecounttype")) %>% tibble()
vstresultsdf <- merge(ddsfinal, vstcountsfinal, by = c("gene_id", "tecounttype")) %>% tibble()

####### GOT HERE


write_tsv(resultsdf, outputs$resultsdf)
write_tsv(vstresultsdf, outputs$vstresultsdf)
