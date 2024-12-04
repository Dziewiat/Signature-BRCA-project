# Install deseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Load deseq2 library
library(DESeq2)

# Choose to filter data
filter <- "ORIGINAL"

if (filter == "FILTERED") {
  count_matrix_filename <- "count_matrix_filtered.tsv"
  normalized_vst_counts_filename <- "normalized_vst_counts_filtered.tsv"
} else if (filter == "ORIGINAL") {
  count_matrix_filename <- "count_matrix.tsv"
  normalized_vst_counts_filename <- "normalized_vst_counts.tsv"
} else if (filter == "EXCLUDED") {
  count_matrix_filename <- "count_matrix_excluded.tsv"
  normalized_vst_counts_filename <- "normalized_vst_counts_excluded.tsv"
}

# Load count data and metadata
countData <- read.csv(paste("counts", count_matrix_filename, sep = "/"),
                      sep = "\t", header = TRUE, row.names = 1)
colData <- read.csv("counts/info_table.tsv", sep = "\t", header = TRUE,
                    row.names = 1)
head(countData[c(1, 2, 3)])
head(colData)

# Create dds
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ study)

# Normalized vst counts
vst_data <- vst(dds, blind = TRUE)
normalized_vst_counts <- assay(vst_data)
head(normalized_vst_counts[c(1, 2, 3)])

write.table(as.data.frame(normalized_vst_counts),
            file = paste("normalized_counts",
                         normalized_vst_counts_filename, sep = "/"),
            sep = "\t")
