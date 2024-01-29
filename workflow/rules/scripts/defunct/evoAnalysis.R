library(seqinr)
library(dplyr)
library(msa)
library(ape)
library(ggmsa)
library(ggtree)
library(cowplot)
library(httpgd)
library(ggplotify)
library(plotly)
library(Biostrings)
library(phangorn)
save.image()
#make bs image incase the aligment was empty
pdf(snakemake@output[["alignment"]], 12, 5)
print(plot(1,1))
dev.off()

try({{
alnFile <- snakemake@input[["alnWithConsensus"]]
aln <- read.dna(alnFile, "fasta")
pdf(snakemake@output[["alignment"]], 12, 5)
par(mar = c(5, 13, 4, 1) + .1)
print(image(aln))
dev.off()
}}, silent = TRUE)

# visualize alignment

# #### Plot MSA
# p <- ggmsa(alnFile, 100, 150, char_width = 0.5, seq_name = TRUE) + geom_seqlogo()
# pdf(snakemake@output[["msa"]], 10,2)
# print(p)
# dev.off()
#### Get distances using ape and then build tree with phangorn
# d <- dist.dna(aln, pairwise.deletion = TRUE)

# try({{treeUPGMA <- upgma(d)
# treeNJ <- nj(d)
# tp <- ggplot(treeUPGMA, aes(x, y)) +
#     geom_tiplab() +
#     geom_tree() +
#     theme_tree() +
#     xlim(0, 0.2) +
#     geom_treescale()
# path = paste0(snakemake@params[["basepath"]], "treeUPGMA.pdf")
# pdf(path, 6, 4)
# print(tp)
# dev.off()
# }})
# datat = tidy_msa(alnFile, 180, 208)
# tp + geom_facet(geom = geom_msa, data = alnFile,  panel = 'msa', color = "Chemistry_AA")

# pdf(snakemake@output[["treeMSA"]], 12, 4)
# print(tp)
# dev.off()
# }}, silent = TRUE)
