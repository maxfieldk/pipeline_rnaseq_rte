library("readr")
library("stringr")
library("dplyr")
library("ggplot2")
library("tibble")
library("RColorBrewer")
library("magrittr")
library("cowplot")
library("eulerr")
library("ggVennDiagram")
library("Gviz")
library("GenomicRanges")
library("rtracklayer")
library("trackViewer")
library("org.Hs.eg.db")
####################################
bams <- snakemake@input[["sortedSTARbams"]]
refseq <- snakemake@params[["refseq"]]
genes <- snakemake@params[["genes"]]
l1hs6kbintactbed <- snakemake@params[["l1hs6kbintactbed"]]
repeats <- snakemake@params[["repeats"]]
peptable <- read.csv(snakemake@params[["peptable"]])
samples <- peptable$sample_name
ideodf <- snakemake@params[["ideogram"]]
ideodf <- "/users/mkelsey/data/ref/genomes/hs1/annotations2/ideogramWithStain.bed"
levels <- snakemake@params[["levels"]]
condition_colors <- snakemake@params[["condition_colors"]]
dertes <- read.delim(snakemake@input[["DETEsbyContrast"]], header = TRUE)
outputdir <- snakemake@params[["outputdir"]]
contrast <- snakemake@wildcards[["contrast"]]
rtekind <- snakemake@wildcards[["rtekind"]]
telocaltype <- snakemake@wildcards[["telocaltype"]]
myrange <- as.integer(snakemake@wildcards[["range"]])
#######################################
# getOption("Gviz.scheme")
# ## [1] "default"
scheme <- getScheme()
scheme$GeneRegionTrack$arrowHeadMaxWidth <- 5
scheme$GeneRegionTrack$arrowHeadWidth <- 5
scheme$GeneRegionTrack$col <- "black"
scheme$GeneRegionTrack$fill <- "#FFD58A"
addScheme(scheme, "myScheme")
options(Gviz.scheme = "myScheme")

# scheme$GeneRegionTrack$shape <- c("smallArrow", "arrow")


# tab <- table(mapper$TE)
# getPos <- function(rte, tab) {
#     if (unname(tab[rte]) > 1) {
#         return(NULL)
#     } else {
#         rterow <- mapper %>% filter(TE == rte)
#         val <- rterow["chromosome.start.stop.strand"]
#         chr <- str_split(val, ":")[[1]][1]
#         coords <- str_split(val, ":")[[1]][2]
#         strand <- str_split(val, ":")[[1]][3]
#         start <- str_split(coords, "-")[[1]][1]
#         stop <- str_split(coords, "-")[[1]][2]

#         return(list(chr = chr, start = as.numeric(start), stop = as.numeric(stop), strand = strand))
#     }
# }

genome <- "hs1"
ideodf <- read.delim(ideodf, header = TRUE)
strack <- SequenceTrack("/users/mkelsey/data/ref/genomes/hs1/hs1.sorted.fa")
# localrtetrack <- AnnotationTrack(start = start, width = stop - start, chromosome = chr, strand = strand, id = rte, genome = genome)
ideogramtrack <- IdeogramTrack(genome = genome, bands = ideodf)
genesTrack <- GeneRegionTrack(genes,
    genome = genome,
    showFeatureId = TRUE,
    showID = TRUE,
    name = "RefSeq",
    stacking = "squish",
    collapseTranscripts = TRUE,
    shape = "arrow"
)
repeatTrack <- GeneRegionTrack(repeats,
    name = "RepeatMasker",
    genome = genome,
    showFeatureId = TRUE,
    showID = TRUE,
    stacking = "squish",
    shape = "arrow",
    fill = "#95c8ff"
)

# generate alignemnt tracks for the relevant samples for this contrast
conditions <- gsub("condition_", "", contrast) %>% str_split_1("_vs_")
samples <- peptable[peptable$condition %in% conditions, ]$sample_name
bamsubset <- c()
for (sample in samples) {
    samplebam <- bams[(grepl(sample, bams))]
    bamsubset <- c(bamsubset, samplebam)
}
alignmentTrackList <- list()
for (i in seq(length(bamsubset))) {
    sample <- str_split_1(basename(bamsubset[i]), "\\.")[1]
    condition <- peptable[peptable$sample_name == sample, ]$condition
    fillcolor <- condition_colors[[condition]]
    assign(
        paste0("alTrack", i),
        AlignmentsTrack(bamsubset[i],
            name = sample,
            isPaired = TRUE,
            genome = "hs1",
            type = "coverage",
            fill = fillcolor
        )
    )
    track <- get(paste0("alTrack", i))
    alignmentTrackList[sample] <- track
}

######
dertedf <- dertes[dertes$telocaltype == telocaltype & dertes$contrast == contrast & dertes$Subfamily == rtekind, ]
for (rte in dertedf$teorgenename) {
    dirname <- file.path(outputdir, telocaltype, contrast, rtekind, rte)
    dir.create(dirname, recursive = TRUE)
    rterow <- dertedf[dertedf$teorgenename == rte]
    direction <- rterow$direction
    genome <- "hs1"
    chr <- rterow$chr
    strand <- rterow$strand
    start <- rterow$start
    stop <- rterow$stop
    length <- rterow$length
    gtrack <- GenomeAxisTrack(name = paste0(length, "bp"))

    tracks <- c(gtrack, alignmentTrackList, repeatTrack, genesTrack)
    lengthtracks <- length(tracks)
    ht <- HighlightTrack(
        trackList = tracks,
        start = start, end = stop,
        chromosome = chr,
        inBackground = TRUE,
        fill = "#d5e2ff",
        col = "transparent"
    )

    ylims <- list(c(0, 10), c(0, 27), c(0, 53), c(0, 105), c(0, 305), c(0, 505))
    for (ylim in ylims) {
        ylim <- unlist(ylim)
        basename <- paste0(
            direction, rte, "_", chr, "_",
            start, "_", stop, "_", "range", myrange,
            "ylim", ylim[2], ".pdf"
        )
        try(
        pdf(file.path(dirname, basename), width = 7, height = 10)
        pl <- plotTracks(
            c(ideogramtrack, ht),
            chromosome = chr,
            from = start,
            to = stop,
            ylim = ylim,
            extend.right = myrange,
            extend.left = myrange,
            sizes = c(0.2, rep(0.3, lengthtracks)),
            showFeatureId = TRUE,
            featureAnnotation = "id",
            shape = "fixedArrow",
            collapseTranscripts = TRUE,
            showId = TRUE,
            col = "black",
            fontcolor = "black",
            fontcolor.group = "black",
            background.panel = NULL,
            background.title = "#dadada",
            col.title = "black",
            col.border.title = "black",
            col.axis = "black",
            col = "black",
            showSampleNames = TRUE
        )
        dev.off(), silent = TRUE
        )
    }
}

#####################

x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
