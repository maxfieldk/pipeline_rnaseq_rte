library(tibble)
library(readr)
library(dplyr)
library(magrittr)
library(stringr)


# fqs <- list.files(recursive = TRUE, pattern = "fastq.gz")
# df <- tibble(X1 = fqs)
# df <- df %>% mutate(leading = str_extract(X1, "^.*/"))
# df <- df %>% mutate(rnanumber = str_extract(X1, "RN[S]*[-]*[0-9][0-9]*"))
# df <- df %>% mutate(rnanumber = str_remove(rnanumber, "S"))
# df <- df %>% mutate(rnanumber = str_remove(rnanumber, "-"))
# df <- df %>% mutate(trailing = str_extract(X1, "_R1|_R2"))
# df <- df %>% mutate(new = paste0("rawdata/", rnanumber, trailing, "_001.fastq.gz"))
# df <- df %>% filter(!str_detect(new, "NA"))

# old <- df$X1
# new <- df$new

# write_csv(df, "renaming.csv")


# df[duplicated(new), ] %>% print(n = 100)
# head(old1, 50)
# for (i in 1:length(old1)) {
#     file.rename(
#         from = old1[i],
#         to = new1[i]
#     )
# }


# metadata %>%
#     mutate(sample_name = paste0("RN", `rna number`)) %>%
#     mutate(R1 = paste0(sampleID, "_R1_001.fastq.gz")) %>%
#     mutate(R2 = paste0("rawdata/", sampleID, "_R2_001.fastq.gz")) %>%
#     rename() %>%
#     relocate(sample_name, condition, R1, R2)


# new
# metadata %>% filter(is.na(R1))

# metadata
metadata <- read_csv("/users/mkelsey/data/denis/samplesheet.csv")
fastqpaths <- list()
metadatanew <- tibble(sample_name = character(), condition = character(), conditionlessdetail = character(), R1 = character(), R2 = character(), R1trimmed = character(), R2trimmed = character())
for (number in metadata$`rna number`) {
    print(number)
    fqs <- list.files(recursive = TRUE, pattern = paste0("RN", number, "([^0-9]).*", ".fastq.gz"))
    print(fqs)
    fastqpaths[[as.character(number)]] <- fqs
    if (length(fqs) == 2) {
        metadatanew <- metadatanew %>% add_row(
            sample_name = paste0("RN", number),
            condition = metadata %>% filter(`rna number` == number) %>% pull(group_detail),
            conditionlessdetail = metadata %>% filter(`rna number` == number) %>% pull(group),
            R1 = fqs[grepl("R1_001", fqs)],
            R2 = fqs[grepl("R2_001", fqs)],
            R1trimmed = sprintf("outs/%s/trimmedReads/%s_1.trimmed.fastq.gz", paste0("RN", number), paste0("RN", number)),
            R2trimmed = sprintf("outs/%s/trimmedReads/%s_2.trimmed.fastq.gz", paste0("RN", number), paste0("RN", number)),
        )
    }
}


# removed from analysis sample 34 and 38, unclear which are their associated fastqs.
write_csv(metadatanew, "conf/private/sample_table.csv")
