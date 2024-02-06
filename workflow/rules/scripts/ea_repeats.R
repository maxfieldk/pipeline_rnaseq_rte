source("~/data/common/myDefaults.r")
library(magrittr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(stringr)
library(cowplot)
library(pathview)
library(ggplot2)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(ggbeeswarm)
library(paletteer)
library(forcats)
library(ggstance)
library(enrichplot)
library(circlize)

# analysis parameters
{
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
                "counttypes" = conf$counttypes,
                "sample_table" = conf$sample_table,
                "contrasts" = conf$contrasts,
                "inputdir" = "results/agg/deseq2/telocal_multi",
                "outputdir" = "results/agg/enrichment_analysis_repeats/telocal_multi",
                "tecounttypes" = c("telescope_multi"),
                "r_annotation_fragmentsjoined" = conf$r_annotation_fragmentsjoined,
                "r_repeatmasker_annotation" = conf$r_repeatmasker_annotation
            ), env = globalenv())
            assign("inputs", list(
                "resultsdf" = "results/agg/repeatanalysis_telescope/resultsdf.tsv"
            ), env = globalenv())
            assign("outputs", list(outfile = "results/agg/enrichment_analysis_repeats/telocal_multi/outfile.txt"), env = globalenv())
        }
    )

    sample_table <- read_csv(params[["sample_table"]])
}


## Load Data and add annotations
resultsdf1 <- read_delim(inputs$resultsdf, delim = "\t")
r_annotation_fragmentsjoined <- read_csv(params$r_annotation_fragmentsjoined)
r_repeatmasker_annotation <- read_csv(params$r_repeatmasker_annotation)
resultsdf <- resultsdf1 %>%
    left_join(r_annotation_fragmentsjoined) %>%
    left_join(r_repeatmasker_annotation)


ontologies <- colnames(r_repeatmasker_annotation)[!(colnames(r_repeatmasker_annotation) %in% c("gene_id", "family"))]
small_ontologies <- ontologies[grepl("subfamily", ontologies)]
big_ontologies <- ontologies[!grepl("subfamily", ontologies)]



gse_results <- list()
core_enrichments_for_plot <- list()

for (contrast in params[["contrasts"]]) {
    contrast_of_interest <- contrast
    contrast_level_1 <- contrast_of_interest %>%
        str_split("_") %>%
        unlist() %>%
        .[4]
    contrast_level_2 <- contrast_of_interest %>%
        str_split("_") %>%
        unlist() %>%
        .[2]
    contrast_stat <- paste0("stat_", contrast_of_interest)
    contrast_padj <- paste0("padj_", contrast_of_interest)
    contrast_log2FoldChange <- paste0("log2FoldChange_", contrast_of_interest)

    # PREP RESULTS FOR GSEA
    res <- resultsdf %>%
        filter(tecounttype == tecounttype) %>%
        arrange(-!!sym(contrast_stat))
    ordered_by_stat <- setNames(res %>% pull(!!sym(contrast_padj)), res$gene_id) %>% na.omit()

    for (ontology in ontologies) {
        tryCatch(
            {
                genesets <- resultsdf %>%
                    select(!!sym(ontology), gene_id) %>%
                    filter(!!sym(ontology) != "Other")
                gse <- GSEA(ordered_by_stat, TERM2GENE = genesets, eps = 0, maxGSSize = 10000, minGSSize = 1)
                gse_results[[contrast]][[ontology]] <- gse

                genesettheme <- theme_gray() + theme(axis.text.y = element_text(colour = "black"))
                p <- dotplot(gse, showCategory = 20) + ggtitle(paste("GSEA", contrast, sep = " ")) + genesettheme + mytheme
                mysave(sprintf("%s/%s/gsea/%s/dotplot.png", params[["outputdir"]], contrast, ontology), w = 4, h = 6, res = 300)

                p <- ridgeplot(gse) + ggtitle(paste("GSEA", contrast, sep = " ")) + genesettheme + mytheme
                mysave(sprintf("%s/%s/gsea/%s/ridgeplot.png", params[["outputdir"]], contrast, ontology), w = 4, h = 6, res = 300)


                tryCatch(
                    {
                        for (num in c(5, 10, 15, 30)) {
                            df <- arrange(gse, -abs(NES)) %>%
                                group_by(sign(NES)) %>%
                                slice(1:num)
                            df <- df@result

                            p <- ggplot(df, aes(NES, fct_reorder(Description, NES), fill = p.adjust)) +
                                geom_col(orientation = "y") +
                                scale_fill_continuous(low = "red", high = "blue", guide = guide_colorbar(reverse = TRUE)) +
                                mytheme +
                                theme(axis.text.y = element_text(colour = "black")) +
                                ylab(NULL)

                            mysave(sprintf("%s/%s/gsea/%s/nes%s.png", params[["outputdir"]], contrast, ontology, num), w = 4, h = min(num, 7), res = 300)
                        }
                    },
                    error = function(e) {
                        print("")
                    }
                )

                tryCatch(
                    {
                        p <- gseaplot2(gse, geneSetID = "SAUL_SEN_MAYO", pvalue_table = TRUE, subplots = 1:2, ES_geom = "line")
                        mysave(sprintf("%s/%s/gsea/%s/senmayo_gsea.png", params[["outputdir"]], contrast, ontology), w = 8, h = 3, res = 300)
                    },
                    error = function(e) {
                        print("")
                    }
                )
            },
            error = function(e) {
                print(e)
            }
        )
    }
}

# GSEA


# annotationForR <- read_csv("/users/mkelsey/data/ref/genomes/hs1/annotations4/hs1.repeatMasker.gtf.fragmentationaccountedfor.rformatted.csv")
# annosmall <- annotationForR %>% sample_n(1000)

# gene_id <- annotationForR$gene_id
# str_count(gene_id, "/") %>% table()
# gene_id[str_count(gene_id, "/") == 4] %>% head()
# gene_id[str_count(gene_id, "/") == 3] %>% head()
# gene_id[str_count(gene_id, "/") == 2] %>% head()

# gene_id[str_count(gene_id, "/") == 2] %>% head()
# annotationForR %>%
#     filter(grepl("L1HS", gene_id)) %>%
#     select(family, gene_id)
# annotationForR %>%
#     tail() %>%
#     select(family, gene_id)

# rtes <- annotationForR %>% filter(grepl("LINE|SINE|ERV", family))

# gene_id <- rtes %>% pull(gene_id)
# str_count(gene_id, "/") %>% table()



# ontology1df <- rtes %>%
#     mutate(term = str_split(gene_id, "/") %>%
#         map_chr(1)) %>%
#     select(term, gene_id) %>%
#     rename(gene = gene_id)

# ontology2df <- rtes %>%
#     mutate(term = str_split(gene_id, "/") %>%
#         map_chr(2)) %>%
#     select(term, gene_id) %>%
#     rename(gene = gene_id)

# ontology3df <- rtes %>%
#     mutate(term = str_split(gene_id, "/") %>%
#         map_chr(3)) %>%
#     select(term, gene_id) %>%
#     rename(gene = gene_id)

# ontologies <- list(
#     ontology1df, ontology2df, ontology3df
# )

# ontology1df <- rtes %>%
#     mutate(term = str_split(gene_id, ":") %>%
#         map_chr(4)) %>%
#     select(term, gene_id) %>%
#     rename(gene = gene_id)

# ontology2df <- rtes %>%
#     mutate(term = str_split(gene_id, ":") %>%
#         map_chr(3)) %>%
#     select(term, gene_id) %>%
#     rename(gene = gene_id)

# ontology3df <- rtes %>%
#     mutate(term = str_split(gene_id, ":") %>%
#         map_chr(2)) %>%
#     select(term, gene_id) %>%
#     rename(gene = gene_id)



# for (ontology in small_ontologies) {
#     ontology_groups <- r_annotation_families %>%
#         pull(!!sym(ontology)) %>%
#         unique()
#     ontology_groups <- ontology_groups[ontology_groups != "Other"]
#     for (group in ontology_groups) {

#     }

# }
