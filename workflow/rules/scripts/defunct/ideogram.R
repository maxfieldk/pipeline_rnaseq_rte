library("readr")
library("stringr")
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
library("RIdeogram")
library("dplyr")

save.image()

log <- file(snakemake@log[[1]], open = "wt")
sink(log)

####### functions
get_ideogram_labels <- function(
    rtedf,
    elementlist,
    directions = c("UP", "DOWN"),
    contrast = "condition_LSEN_vs_PRO",
    tecountstype = "telocal_multi",
    namedmarkerlist = c("L1" = "circle", "Alu" = "box", "ERVK" = "triange", "L1HS" = "circle", "AluY" = "box", "HERVK-int" = "triangle"),
    namedcolorlist = c(
        "L1UP" = "8c00ff",
        "AluUP" = "ff9d00",
        "ERVKUP" = "ff7be2",
        "L1HSUP" = "8c00ff",
        "AluYUP" = "ff9d00",
        "HERVK-intUP" = "ff7be2",
        "L1DOWN" = "00ff62",
        "AluDOWN" = "0be7e4",
        "ERVKDOWN" = "8c00ff",
        "L1HSDOWN" = "00ff62",
        "AluYDOWN" = "0be7e4",
        "HERVK-intDOWN" = "8c00ff"
    )) {
    tempdf <- rtedf %>% dplyr::select(c(Subfamily, chr, start, stop, direction))
    colnames(tempdf) <- c("Kind", "Chr", "Start", "End", "Direction")
    tempdf <- tempdf %>%
        filter(Kind %in% elementlist) %>%
        filter(Direction %in% directions)
    tempdf <- tempdf %>%
        mutate(Shape = namedmarkerlist[Kind]) %>%
        mutate(color = unname(unlist(namedcolorlist[paste0(Kind, Direction)]))) %>%
        mutate(Type = paste(Kind, Direction, sep = " "))
    # now the df is in a format that RIdeogram can read
    df <- tempdf[, c("Type", "Shape", "Chr", "Start", "End", "color")]
    return(df)
}

makeIdeogram <- function(
    #' Be sure to have svgsavepath end kn .svg
    rtedf,
    elementlist,
    karyotype,
    genedensity,
    directions = c("UP", "DOWN"),
    svgsavepath = "chromosome.svg",
    contrast = "condition_LSEN_vs_PRO", tecountstype = "telocal_multi",
    namedmarkerlist = c("L1" = "circle", "Alu" = "box", "ERVK" = "triange", "L1HS" = "circle", "AluY" = "box", "HERVK-int" = "triangle"),
    namedcolorlist = c("L1" = "ff0000", "Alu" = "ffa500", "ERVK" = "0000ff", "L1HS" = "ff0000", "AluY" = "ffa500", "HERVK-int" = "0000ff")) {
    karyotypedf <- read.delim(karyotype, sep = "\t")

    genedensity <- read.delim(genedensity, header = FALSE)
    colnames(genedensity) <- c("Chr", "Start", "End", "Value")
    genedensity$Start <- as.numeric(genedensity$Start)
    genedensity$End <- as.numeric(genedensity$End)
    genedensity$Value <- as.integer(genedensity$Value)

    dertes <- read.delim(rtedf, sep = "\t", header = TRUE)
    rownames(dertes) <- 1:nrow(dertes)

    labels <- get_ideogram_labels(dertes, elementlist,
        directions = directions,
        contrast = contrast,
        tecountstype = tecountstype,
        namedcolorlist = namedcolorlist
    )
    ideogram(
        karyotype = karyotypedf,
        overlaid = genedensity,
        label = labels,
        label_type = "marker",
        output = svgsavepath
    )
    convertSVG(svgsavepath, gsub(".svg", ".png", svgsavepath), device = "png")
}

########
rtedf <- snakemake@input[["DETEsbyContrast"]]
genedensitylist <- snakemake@params[["genedensity"]]
rtestoplot <- snakemake@params[["rtestoplot"]]
namedcolorlist <- snakemake@params[["namedcolorlist"]]
namedmarkerlist <- snakemake@params[["namedmarkerlist"]]
karyotype <- snakemake@params[["karyotype"]]
outputdir <- snakemake@params[["outputdir"]]
telocaltypes <- snakemake@params[["telocaltypes"]]
contrasts <- snakemake@params[["contrasts"]]

for (telocaltype in telocaltypes) {
    for (contrast in contrasts) {
        for (rte in rtestoplot) {
            savepath <- paste(outputdir, telocaltype, contrast, paste0(rte, "_ideogram.svg"), sep = "/")
            genedensity <- genedensitylist[[rte]]
            elementlist <- c(rte)
            makeIdeogram(rtedf, elementlist, karyotype, genedensity = genedensity, namedcolorlist = namedcolorlist, svgsavepath = savepath)
        }
    }
}

x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
