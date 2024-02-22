library("DESeq2")
library("readr")
library("pheatmap")
library("ggplot2")
library(GenomicRanges)
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
library(ggplot2)
library("readr")
library("magrittr")
library(rstatix)
library(ggpubr)


### functions
plotSave <- function(path, plot, width = 6, height = 6) {
    dir.create(dirname(path), recursive = TRUE)
    png(path, width = width, height = height, units = "in", dpi = 300)
    print(plot)
    dev.off()
}

sample_table <- read_csv("conf/private/sample_table.csv")
samples <- sample_table %$% sample_name

conditions <- unique(sample_table$condition)
condition1samples <- sample_table[sample_table$condition == conditions[1], ]$sample_name
condition2samples <- sample_table[sample_table$condition == conditions[2], ]$sample_name


{#plotting params
my_palette <- paletteer::paletteer_d("ggsci::nrc_npg")

bluelight <- "#8399c9"
blue <- "#3C5488FF"
bluedark <- "#273759"
orangelight <- "#F39B7FFF"
orange <- "#e94716"
orangedark <- "#75240b"
mycolor <- "#00A087FF"
mythemesamples <- list(
    scale_fill_manual(values = rep(my_palette[4:5], each = 3)),
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
)
mythemesamplesdif <- list(
    scale_fill_manual(values = c(bluelight, blue, bluedark, orangelight, orange, orangedark)),
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
)
mythemeconditions <- list(
    scale_color_manual(values = my_palette[5:4]),
    scale_fill_manual(values = my_palette[5:4]),
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
)
mythemecontrast <- list(
    scale_fill_manual(values = my_palette[c(1, 2, 9)]),
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
)
mythemecontrastrev <- list(
    scale_fill_manual(values = my_palette[c(2,1, 9)]),
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
)
}



firstdf <- read_delim(paste0("outs/", samples[1], "/", samples[1], "rtesandgenes.counts.txt"), delim = "\t")
geneid <- firstdf %$% Geneid
bounddf <- data.frame(Geneid = firstdf[, 1])
for (sample in samples) {
    bounddf <- cbind(bounddf, read.delim(paste0("outs/", sample, "/", sample, "rtesandgenes.counts.txt"))[, 2])
}

colnames(bounddf) <- c("Geneid", samples)
bounddf <- bounddf %>% tibble()






### inputs

cts <- dplyr::select(bounddf, -Geneid)
cnames <- colnames(cts)

coldata <- sample_table

contrasts <- "condition_DMD_vs_Ctrl"
contrast <- "condition_DMD_vs_Ctrl"

levels <- c("Ctrl", "DMD")

outputdir <- "results/agg/deseq2/featurecounts"
dir.create(outputdir)
###
condition <- coldata$condition
dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~condition
)
# used to filter here
# this sets prol as the reference level since its first in the vector
dds$condition <- factor(dds$condition, levels = levels)

####
dds <- DESeq(dds)
save.image()
keep <- rowSums(counts(dds)) >= 10
ddsfiltered <- dds[keep, ]
####
resultsNames(dds) # lists the coefficients

####
# for pca ill use the ddsfiltered but will keep the unfiltered for repeat analysis purposes
# most repeats just have a couple counts - these should not be thrown away, they can be aggregated
length(geneid)
length(dds)
maskrte = !grepl("/", geneid)
ddsnoRTEs = dds[maskrte,]
keep <- rowSums(counts(ddsnoRTEs)) >= 10
ddsnoRTEsfiltered = ddsnoRTEs[keep, ]
vst <- assay(vst(ddsnoRTEsfiltered))
vst_full <- vst(ddsnoRTEsfiltered)
sampleDists <- dist(t(vst))
counttablesizenormed <- counts(dds, normalized = T)

## PCA plots
pcaObj <- pca(vst, metadata = colData(dds), removeVar = 0.1)

screep <- screeplot(pcaObj, title = "") +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

dir.create(paste(outputdir, "plots", sep = "/"))
png(paste(outputdir, "plots", "screeplot.png", sep = "/"), width = 5, height = 6, units = "in", res = 300)
print(screep)
dev.off()

loadingsp <- plotloadings(pcaObj,
    components = getComponents(pcaObj, seq_len(3)),
    rangeRetain = 0.045, labSize = 4
) +
    theme(legend.position = "none") +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, "plots", "loadings.png", sep = "/"), width = 6, height = 6, units = "in", res = 300)
print(loadingsp)
dev.off()

pcaplot <- biplot(pcaObj,
    showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
)
pcap <- pcaplot +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, "plots", "pcaplot.png", sep = "/"), width = 7, height = 6, units = "in", res = 300)
print(pcap)
dev.off()

pcaplot34 <- biplot(pcaObj,
    x = "PC3", y = "PC4", showLoadings = FALSE, gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0,
    colby = "condition", legendPosition = "right",
    labSize = 5, pointSize = 5, sizeLoadingsNames = 5
)
pcap34 <- pcaplot34 +
    theme_cowplot() +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))

png(paste(outputdir, "plots", "pcaplotPC34.png", sep = "/"), width = 7, height = 6, units = "in", res = 300)
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

png(paste(outputdir, "plots", "PCAgrid.png", sep = "/"), width = 16, height = 6, units = "in", res = 300)
print(p)
dev.off()

# png("eigencorr.png", width=10, height=8)
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
png(paste(outputdir, "plots", "pcaplot_statelipse.png", sep = "/"), width = 8, height = 8, units = "in", res = 300)
print(pcaplot_statelipse)
dev.off()


############

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst_full$condition, vst_full$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(paste(outputdir, "plots", "heatmapplot.png", sep = "/"), width = 10, height = 8, units = "in", res = 300)
pheatmap::pheatmap(sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors
)
dev.off()
####

res <- results(dds, name = contrast)
# or to shrink log fold changes association with condition:
resLFC <- lfcShrink(dds, coef = contrast, type = "apeglm")
# Save the plot as a png file
library(EnhancedVolcano)
png(paste(outputdir, "deplot.png", sep = "/"), width = 10, height = 8, units = "in", res = 300)
z <- EnhancedVolcano(res,
    lab = rownames(res),
    selectLab = c(""),
    title = "Aged vs Young",
    drawConnectors = TRUE,
    x = "log2FoldChange",
    y = "padj",
    ylab = expression(-Log[10] ~ P["adj"]),
    legendPosition = "none",
    pCutoff = 0.05,
    xlim = c(-10, 10),
    ylim = c(-1, 15)
) +
    theme(axis.line = element_blank(), aspect.ratio = 1, panel.border = element_rect(color = "black", linetype = 1, linewidth = 1, fill = NA))
print(z)
dev.off()

resOrdered <- res[order(res$pvalue), ]
resSig <- subset(resOrdered, padj < 0.1)
counttablesizenormed <- counts(dds, normalized = T)
rld <- rlog(dds, blind = FALSE)

active_conditions <- str_split_1(gsub("condition_", "", contrast), "_vs_")
tmp <- coldata %>% filter(str_detect(condition, active_conditions[1]) | str_detect(condition, active_conditions[2]))
df <- data.frame(condition = tmp$condition)
rownames(df) <- tmp$sample_name
topexpressed <- order(rowMeans(counttablesizenormed[, rownames(df)]),
    decreasing = TRUE
)[1:40]
topSig <- resOrdered %>%
    head(n = 40) %>%
    rownames()
topVarGenes <- head(order(rowVars(assay(rld)[, rownames(df)]), decreasing = TRUE), 500)
topVarGenes40 <- head(order(rowVars(assay(rld)[, rownames(df)]), decreasing = TRUE), 40)
interestingsets <- list("Top 40 Variable Genes" = topVarGenes40, "Top Highly Expressed Genes" = topexpressed, "Top Differentially Expressed Genes" = topSig, "Top Variable Genes" = topVarGenes)


for (name in names(interestingsets)) {
    for (scale in c("none", "row")) {
        for (binary in c(TRUE, FALSE)) {
            if (binary) {
                blabel <- "rownames"
            } else {
                blabel <- "norownames"
            }
            set <- interestingsets[[name]]
            dirname <- file.path(outputdir, counttype, contrast)
            dir.create(dirname, recursive = TRUE)
            filename <- paste0(name, blabel, scale, ".png")
            path <- file.path(dirname, filename)
            png(path, width = 10, height = 8, units = "in", res = 300)
            pheatmap(assay(rld)[set, rownames(df)], annotation_col = df, scale = scale, show_rownames = binary, main = name)
            dev.off()
        }
    }
}


write.csv(as.data.frame(resOrdered), file = paste(outputdir, "results.csv", sep = "/"))
write.csv(as.data.frame(resSig), file = paste(outputdir, "resultsSig.csv", sep = "/"))
write.csv(as.data.frame(counttablesizenormed), file = paste(outputdir, "counttablesizenormed.csv", sep = "/"))
write.csv(as.data.frame(assay(rld)), file = paste(outputdir, "rlogcounts.csv", sep = "/"))
}

counttablesizenormed <- counts(dds, normalized = T)
rld <- rlog(dds, blind = FALSE)

tmp <- coldata
df <- data.frame(condition = tmp$condition)
rownames(df) <- tmp$sample_name
topexpressed <- order(rowMeans(counttablesizenormed),
    decreasing = TRUE
)[1:40]
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 500)
topVarGenes40 <- head(order(rowVars(assay(rld)), decreasing = TRUE), 40)
interestingsets <- list("Top 40 Variable Genes" = topVarGenes40, "Top Highly Expressed Genes" = topexpressed, "Top Variable Genes" = topVarGenes)


for (name in names(interestingsets)) {
    for (scale in c("none", "row")) {
        for (binary in c(TRUE, FALSE)) {
            if (binary) {
                blabel <- "rownames"
            } else {
                blabel <- "norownames"
            }
            set <- interestingsets[[name]]
            dirname <- file.path(outputdir, counttype, "plots")
            dir.create(dirname, recursive = TRUE)
            filename <- paste0(name, blabel, scale, ".png")
            path <- file.path(dirname, filename)
            png(path, width = 10, height = 8, units = "in", res = 300)
            pheatmap(assay(rld)[set, ], annotation_col = df, scale = scale, show_rownames = binary, main = name)
            dev.off()
        }
    }
}
}

x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)







#######
df <- counttablesizenormed %>% as.data.frame %>% tibble()
df$geneid <- geneid
padj <- res@listData$padj
df$padj <- padj
rmfragments <- read_csv("/users/mkelsey/data/ref/genomes/hs1/annotations4/hs1.repeatMasker.gtf.fragmentationaccountedfor.rformatted.csv")
j = left_join(df, rmfragments, join_by(geneid == gene_id))
tdf = j %>% pivot_longer(-c(colnames(j)[!(colnames(j) %in% samples)]), names_to = "sample", values_to = "counts")
tdf = tdf %>% mutate(condition = ifelse(sample %in% condition1samples, conditions[1], conditions[2]))
write_csv(tdf, "results/agg/tdf.csv")


# # rmfragments = rm %>% group_by(gene_id) %>% summarise(chr = first(chr), source = first(source), feature = first(feature),  start = min(start), end = max(end), strand = first(strand), frame = first(frame), family = first(family), element_start = min(element_start), element_end = max(element_end), pctdiv = sum(pctdiv*length)/sum(length), length = sum(length), num_fragments = n())
# write_csv(rmfragments, "/users/mkelsey/data/ref/genomes/hs1/annotations4/hs1.repeatMasker.gtf.fragmentationaccountedfor.rformatted.csv")
#annotate intact
# rml1hs = rmfragments %>% filter(grepl("L1HS", gene_id))
# rml1hsgrs <- GRanges(
#     seqnames = rml1hs$chr,
#     ranges = IRanges(start = rml1hs$start, end = rml1hs$end),
#     strand = rml1hs$strand,
#     gene_id = rml1hs$gene_id
# )


# l1hsintactpath <- "/oscar/data/jsedivy/mkelsey/ref/genomes/hs1/annotations3/RTE/l1hsintact.bed"
# l1hsintactdf <- read_delim(l1hsintactpath, col_names = FALSE)
# l1hsintact <- GRanges(
#     seqnames = l1hsintactdf$X1,
#     ranges = IRanges(start = l1hsintactdf$X2, end = l1hsintactdf$X3),
#     strand = l1hsintactdf$X6,
#     name = l1hsintactdf$X4,
#     uid = paste(l1hsintactdf$X1, l1hsintactdf$X2, l1hsintactdf$X3, l1hsintactdf$X6, sep = "_"),
#     region = l1hsintactdf$X11
# )

# sbo = subsetByOverlaps(rml1hsgrs, l1hsintact, minoverlap = 5000, ignore.strand = TRUE)
# rmfragments$intact = "n"
# rmfragments[rmfragments$gene_id %in% sbo$gene_id, "intact"] = "y"


l1hs <- tdf %>% filter(grepl("L1HS", geneid)) %>% filter(num_fragments ==1)
AluY <- tdf %>% filter(grepl("AluY", geneid))%>% filter(num_fragments ==1)
HERVK <- tdf %>% filter(grepl("HERVK", geneid))%>% filter(num_fragments ==1)
l1 <- tdf %>% filter(grepl("L1", geneid))
alu <- tdf %>% filter(grepl("Alu", geneid))
herv <- tdf %>% filter(grepl("HERV", geneid))

pvp <- function(df, path="temp.png", title="", width=5, height=4) {
p <- df %>% group_by(geneid, condition) %>% summarise(mean(counts), padj = first(padj)) %>% pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
    ggplot() +
    geom_point(aes(x = Ctrl, y = DMD, color = padj<0.05), alpha = 0.4) +
    scale_color_manual(values = c("black", "red")) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Ctrl Counts", y = "DMD Counts") +
    ggtitle(title) +
    panel_border(color = "black", linetype = 1, remove = FALSE)
png(path, width = width, height = height, units = "in", res = 300)
print(p)
dev.off()
return(p)
}

pvp(l1hs %>% filter(length > 6000) %>% filter(intact== "y"), "results/agg/plots/l1hsintact.png", "L1HS (>6kb) Intact")
pvp(l1hs %>% filter(length > 6000), "results/agg/plots/l1hs.png", "L1HS (>6kb) Counts")
pvp(AluY %>% filter(length > 300), "results/agg/plots/aluY.png", "AluY (>300bp) Counts")
pvp(HERVK %>% filter(length > 5000), "results/agg/plots/HERVK.png", "HERVK (>5kb) Counts")


stripp <- function(df, path="temp.png",title="", width=5, height=4) {
stat.test = df %>% group_by(sample) %>% summarise(sum(counts), condition = first(condition)) %>%
    t_test(`sum(counts)` ~ condition) %>%
    add_significance() %>% add_xy_position(x = "condition")
stat.test$p.format <- p_format(
  stat.test$p, accuracy = 0.001,
  leading.zero = TRUE
  )

p <- df %>% group_by(sample) %>% summarise(sum(counts), condition = first(condition)) %>%
    ggplot(aes(x = condition, y = `sum(counts)`)) +
    stat_summary(fun = mean, geom = "bar", aes(fill = condition)) +
    geom_point() +
    panel_border(color = "black", linetype = 1, remove = FALSE) +
    labs(x = "", y = "Counts") +
    ggtitle(title) +
    mythemecontrastrev +
    stat_pvalue_manual(stat.test, label = "p.format", bracket.nudge.y = 10) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

png(path, width = width, height = height, units = "in", res = 300)
print(p)
dev.off()
return(p)
}

stripp(l1hs %>% filter(length > 6000) %>% filter(intact== "y"), "results/agg/plots/stripplotl1hsintact.png", "L1HS (>6kb) Intact Counts")
stripp(l1hs %>% filter(length > 6000), "results/agg/plots/stripplotl1hs.png", "L1HS (>6kb) Counts")
stripp(AluY %>% filter(length > 250), "results/agg/plots/stripplotaluY.png", "AluY (>300bp) Counts")
stripp(HERVK %>% filter(length > 5000), "results/agg/plots/stripplothervk.png", "HERVK (>5kb) Counts")


stripp(l1, "results/agg/plots/stripplotl1.png", "L1 Counts")
stripp(alu, "results/agg/plots/stripplotalu.png", "Alu Counts")
stripp(herv, "results/agg/plots/stripplotherv.png", "HERV Counts")


dep <- function(df, path="temp.png", title="", width=5, height=4) {
p <- df %>% group_by(geneid, condition) %>%
summarise(mean(counts), padj = first(padj)) %>%
pivot_wider(names_from = condition, values_from = `mean(counts)`) %>%
mutate(direction = ifelse(padj < 0.05, ifelse(DMD > Ctrl, "UP", "DOWN"), "NS")) %>%
mutate(direction = ifelse(is.na(direction), "NS", direction)) %>%
mutate(direction = direction %>% factor(levels = c("UP", "DOWN", "NS"))) %>%
ungroup() %>%
mutate(n = n()) %>%
group_by(direction, n) %>%
summarise(count = n()) %>%
mutate(fraction = count / n) %>%
filter(direction != "NS") %>%
ggplot() + 
ggtitle(title) +
labs(x = "", y = "Fraction Differentially Expressed") +
geom_col(aes(x = direction, fill = direction, y = fraction)) +
mythemecontrastrev
png(path, width = width, height = height, units = "in", res = 300)
print(p)
dev.off()
return(p)
}

dep(l1hs %>% filter(length > 6000) %>% filter(intact== "y"), "results/agg/plots/deplotl1hsintact.png", "L1HS (>6kb) Intact")
dep(l1hs %>% filter(length > 6000), "results/agg/plots/deplotl1hs.png", "L1HS (>6kb)")
dep(AluY %>% filter(length > 250), "results/agg/plots/deplotaluY.png", "AluY (>300bp)")
dep(HERVK %>% filter(length > 5000), "results/agg/plots/deplothervk.png", "HERVK (>5kb)")

