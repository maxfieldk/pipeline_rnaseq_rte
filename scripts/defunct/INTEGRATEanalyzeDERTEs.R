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
library("ggrepel")
library("grid")

telocaltypes <- snakemake@params[["telocaltypes"]]
contraststocompare <- snakemake@params[["contraststocompare"]]
peptable <- read.csv(snakemake@params[["peptable"]])
outputdir <- snakemake@params[["outputdir"]]


##############
makeLegendGrob <- function(labels, fill = c("blue", "red"), cex = 1) {
    return(legendGrob(labels,
        pch = 21,
        gp = gpar(
            col = "black",
            fill = fill,
            cex = cex
        )
    ))
}

euPlot <- function(fit, contrasts, fill = c("blue", "red"), alpha = 0.7, legendcex = 1, legend = TRUE, labels = FALSE, main = "") {
    if (legend == TRUE) {
        legd <- makeLegendGrob(contrasts, fill, legendcex)
        p <- plot(
            fit,
            main = main,
            fills = list(fill = fill, alpha = alpha),
            legend = FALSE,
            quantities = list(col = "black"),
            labels = labels
        )
        row <- plot_grid(p)
        plot <- plot_grid(row, legd, nrow = 2, rel_heights = c(1, 0.25))
    } else {
        plot(
            fit,
            main = main,
            fills = list(fill = fill, alpha = alpha),
            legend = FALSE,
            quantities = list(col = "black"),
            labels = labels
        )
    }
}

###########
datasets <- peptable$batch %>% unique()
dedfholder <- list()
for (dataset in datasets) {
    try(
        dedfholder[[dataset]] <- read.table(
            paste(peptable[peptable$batch == dataset, ][1, "basepath"],
                "results/agg/repeatanalysis/allactiveDETEs.tsv",
                sep = "/"
            )
        ),
        silent = TRUE
    )
}
rtekinds <- dedfholder[[1]]$V1 %>% unique()


for (telocaltype in telocaltypes) {
    for (rtekind in rtekinds) {
        for (direction in c("UP", "DOWN")) {
            filtereddedfholder <- list()
            for (dataset in names(dedfholder)) {
                tempdf <- dedfholder[[dataset]][dedfholder[[dataset]]$V1 == rtekind & dedfholder[[dataset]]$V7 == direction & dedfholder[[dataset]]$V8 %in% contraststocompare & dedfholder[[dataset]]$V9 == telocaltype, ]$V2
                filtereddedfholder[[dataset]] <- tempdf
            }

            main <- paste(rtekind, direction, gsub("condition_", "", contraststocompare[1]))
            fit <- euler(filtereddedfholder, shape = "ellipse")
            p <- euPlot(fit, names(filtereddedfholder), main = main, legend = TRUE, labels = TRUE, legendcex = 1)

            pdf(paste(outputdir, telocaltype, paste0("SHARED", rtekind, direction, ".pdf"), sep = "/"), width = 4, height = 4)
            print(p)
            dev.off()
        }
    }
}



x <- data.frame()
write.table(x, file = snakemake@output[["outfile"]], col.names = FALSE)
