args <- commandArgs(trailingOnly=TRUE)
files <- args

library(GSEABase)
library(GOstats)

dir.create("go_results", showWarnings=FALSE)

for (f in files) {
  res <- read.csv(f, row.names=1)

  up   <- rownames(res[res$log2FoldChange > 0,])
  down <- rownames(res[res$log2FoldChange < 0,])

  ref <- rownames(res)

  gsc <- GeneSetCollection(lapply(ref, function(x)
      GeneSet(setName=x, geneIds=x)))

  go_enriched <- function(gene_ids) {
    params <- GSEAGOHyperGParams(
      geneSetCollection = gsc,
      geneIds = gene_ids,
      universeGeneIds = ref,
      ontology = "BP",
      pvalueCutoff = 0.05,
      conditional = FALSE,
      testDirection = "over"
    )
    summary(hyperGTest(params))
  }

  if (length(up)>0)
    write.csv(go_enriched(up),
              paste0("go_results/", basename(f), "_GO_up.csv"))

  if (length(down)>0)
    write.csv(go_enriched(down),
              paste0("go_results/", basename(f), "_GO_down.csv"))
}

