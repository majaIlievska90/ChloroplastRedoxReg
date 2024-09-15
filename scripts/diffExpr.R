##Differential expression analysis with DESeq2.
library(DESeq2)
library(apeglm)
library(GSEABase)
library(Category)
library(GOstats)

# Load and filter data for 1h and 4w
load("dds_forDeseq_transcrSummarized.RData")
dds1H <- dds[rowSums(counts(dds)) >= 10, ]
rm(dds)

load("4_weeks/dds_forDeseq.RData")
dds4w <- dds[rowSums(counts(dds)) >= 10, ]
rm(dds)

# Prepare counts matrix and colData for 4w
counts.mat <- counts(dds4w)[, c(1,4:10,2,11:18,3,19:21)]
colnames(counts.mat)[1:3] <- paste("merged", colnames(counts.mat)[1:3], sep = "-")
coldata <- data.frame(condition = factor(dds4w$condition))
#dds4w2 <- DESeqDataSetFromMatrix(counts.mat, colData = coldata, design = ~ condition)

# Choose dataset: 1h or 4w
dds_dif <- dds4w  

# Differential Expression Analysis
diff_expr <- function(dds, ref_cond, test_cond) {
  dds$condition <- relevel(dds$condition, ref = ref_cond)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", test_cond, ref_cond))
  res_shrink <- lfcShrink(dds, coef = paste0("condition_", test_cond, "_vs_", ref_cond), type = "apeglm")
  subset(res_shrink, padj < 0.05)
}

# 470nm vs 690nm
resSig_470_690 <- diff_expr(dds_dif, "470nm", "690nm")

# 660nm vs 690nm
resSig_660_690 <- diff_expr(dds_dif, "660nm", "690nm")

# Group comparison for different conditions
group_comparison <- function(dds, group_labels, g1, g2) {
  coldata <- cbind(colData(dds), comparison = factor(group_labels))
  dds <- DESeqDataSetFromMatrix(counts(dds), colData = coldata, design = ~ comparison)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("comparison", g2, g1))
  res_shrink <- lfcShrink(dds, coef = paste0("comparison_", g2, "_vs_", g1), type = "apeglm")
  subset(res_shrink, padj < 0.05)
}

# Comparisons
# 420+470 vs 520, 560, 630, 660
resSig_group1 <- group_comparison(dds_dif, c(rep("g1", 6), rep("g2", 12), rep("none", 3)), "g1", "g2")

# 520+560 vs rest
resSig_group2 <- group_comparison(dds_dif, c(rep("g1", 6), rep("g2", 6), rep("g1", 9)), "g1", "g2")

# 660+630 vs others
resSig_group3 <- run_group_comparison(dds1H, c(rep("g1", 12), rep("g2", 6), rep("g1", 3), rep("none", 4)), "g1", "g2")


##GO enrichment
# Define reference genes and gene set collection
ref_genes <- rownames(assays(dds1H)$count)
gsc <- geneSetCol(rownames(dds1H))  

# Create lists for results, up and downregulated genes
GSEAGOup <- vector('list', length(resSig))   
GSEAGOdown <- vector('list', length(resSig)) 

# Function for GO enrichment analysis
GO_enrichment <- function(gene_ids, gene_set_collection, ref_genes) {
  params <- GSEAGOHyperGParams(
    name = 'Arabidopsis thaliana GO',
    geneSetCollection = gene_set_collection,
    geneIds = gene_ids,
    universeGeneIds = ref_genes,
    ontology = 'BP',
    pvalueCutoff = 0.05,
    conditional = FALSE,
    testDirection = 'over'
  )
  summary(hyperGTest(params))
}

# Loop over resSig to perform GO enrichment analysis for up and downmregulated genes
for (i in seq_along(resSig)) {
  up_genes <- rownames(resSig[[i]][resSig[[i]]$log2FoldChange > 0, ])
  down_genes <- rownames(resSig[[i]][resSig[[i]]$log2FoldChange < 0, ])
  
  # GO enrichment for upregulated genes
  if (length(up_genes) > 0) {
    GSEAGOup[[i]] <- GO_enrichment(up_genes, gsc, ref_genes)
  }
  
  # GO enrichment for downregulated genes
  if (length(down_genes) > 0) {
    GSEAGOdown[[i]] <- GO_enrichment(down_genes, gsc, ref_genes)
  }
}