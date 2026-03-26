args <- commandArgs(trailingOnly = TRUE)
file1h <- args[1]
file4w <- args[2]

library(DESeq2)

load(file1h)
dds1h <- dds[rowSums(counts(dds)) >= 10, ]
saveRDS(dds1h, file="dds1h.rds")

load(file4w)
dds4w <- dds[rowSums(counts(dds)) >= 10, ]
saveRDS(dds4w, file="dds4w.rds")
