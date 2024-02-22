log <- file(snakemake@log[[1]], open="wt")
sink(log)
save.image()
telocaltypes = snakemake@params[["telocaltypes"]]

for (telocaltype in telocaltypes) {
    df = read.delim(snakemake@input[[telocaltype]][1])
    coldata = read.csv(snakemake@params[["sample_table"]])

    bounddf = data.frame(Geneid = df[,1])
    for (sample in snakemake@input[[telocaltype]]) {
        bounddf = cbind(bounddf, read.delim(sample)[,2])
    }
    colnames(bounddf) = c("Geneid", coldata$sample_name)
    write.table(bounddf, file = snakemake@output[[telocaltype]], row.names = FALSE, sep = '\t')

}

