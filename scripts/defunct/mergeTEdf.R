df <- read.delim(snakemake@input[["telocal"]][1])
coldata <- read.csv(snakemake@params[["sample_table"]])

bounddf <- data.frame(Geneid = df[, 1])
for (sample in snakemake@input[["telocal"]]) {
    bounddf <- cbind(bounddf, read.delim(sample)[, 2])
}
colnames(bounddf) <- c("Geneid", coldata$sample_name)
write.table(bounddf, file = snakemake@output[["aggcounts"]], row.names = FALSE, sep = "\t")

#################################

# coldata <- read_csv("/users/mkelsey/data/denis/conf/private/sample_table.csv")

# samplepaths <- paste0("/users/mkelsey/data/denis/outs/", coldata$sample_name, "/TElocal/", coldata$sample_name, "_uniq.cntTable")
# df <- read.delim(samplepaths[1])


# bounddf <- data.frame(Geneid = df[, 1])
# for (sample in samplepaths) {
#     bounddf <- cbind(bounddf, read.delim(sample)[, 2])
# }
# colnames(bounddf) <- c("Geneid", coldata$sample_name)
# write.table(bounddf, file = "outs/agg/TElocalCounts_uniq.txt", row.names = FALSE, sep = "\t")
