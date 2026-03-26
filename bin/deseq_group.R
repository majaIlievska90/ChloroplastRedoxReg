library(DESeq2)
library(apeglm)

dds <- readRDS(commandArgs(trailingOnly=TRUE)[1])

comparisons <- list(
  c("470nm","690nm"),
  c("660nm","690nm")
)

dir.create("pairwise_results", showWarnings=FALSE)

for (c in comparisons) {
  ref  <- c[1]
  test <- c[2]

  dds$condition <- relevel(dds$condition, ref=ref)
  dds <- DESeq(dds)

  res <- results(dds, contrast=c("condition", test, ref))
  shr <- lfcShrink(dds, coef=paste0("condition_",test,"_vs_",ref), type="apeglm")

  sig <- subset(shr, padj < 0.05)
  write.csv(as.data.frame(sig),
            file = paste0("pairwise_results/", test, "_vs_", ref, ".csv"))
}
