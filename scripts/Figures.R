############### FIGURE 2B ###############################
#####DOT PLOT OF GO ENRICHMENT###############
##read the enrichment results for all the comparisons; filter 
##


#navigate to the folder 
##setwd
# Load necessary libraries
library(openxlsx)
library(ggplot2)
library(viridis)
library(reshape2)
library(gridExtra)
library(grid)
library(dplyr)
library(ggVennDiagram)
library(VennDiagram)
library(openxlsx)
library(venneuler)
library(dplyr)

# Read data files
core <- read.xlsx("PSI_control_vs_PSII_control.xlsx", cols = 1)$Gene
ox1h <- read.xlsx("PSI60_vs_PSIIc_enrich.xlsx")
ox4m <- read.xlsx("PSI4_vs_PSIIc_enrich.xlsx")
ox30m <- read.xlsx("PSI30_vs_PSIIc_enrich.xlsx")
psi30_psi4 <- read.xlsx("diffExpressedPSI30_vs_PSI4.xlsx")
psi60_psi4 <- read.xlsx("diffExpressedPSI60_vs_PSI4.xlsx")
psii30_psii4 <- read.xlsx("diffExpressedPSII30_vs_PSII4.xlsx")
psii60_psii4 <- read.xlsx("diffExpressedPSII60_vs_PSII4.xlsx")
psi30_psii30 <- read.xlsx("Diff_expressedPSI30vsPSII30.xlsx")
c660_690 <- read.xlsx("condition_660nm_vs_690nm_enrich.xlsx")
ox_red_contr <- read.xlsx("PSI_PSII_c_enrich.xlsx")
red4m <- read.xlsx("PSII4_vs_PSIc_enrich.xlsx")
c660_630 <- read.xlsx("condition_660nm_vs_630nm_enrich.xlsx")
c470_420 <- read.xlsx("condition_470nm_vs_420nm_enrich.xlsx")
c_red_ox <- read.xlsx("condition_reducing1h_vs_oxidizing1h_enrich.xlsx")
c560_520 <- read.xlsx("condition_560nm_vs_520nm_enrich.xlsx")
red30m <- read.xlsx("PSII30_vs_PSIc_enrich.xlsx")
red1h <- read.xlsx("PSII60_vs_PSIc_enrich.xlsx")
c420_470_660 <- read.xlsx("c_red420_470_660_vs_rest.xlsx")

# Create a list of comparisons and their names
comparisons <- list(ox1h, ox30m, ox4m, red1h, red30m, red4m, ox_red_contr, c_red_ox, c660_630, c660_690, c560_520, c470_420)
names <- list("ox_1h", "ox_30m", "ox_4m", "red_1h", "red_30m", "red_4m", "ox_red_contr", "c_red_ox", "c660_630", "c660_690", "c560_520", "c470_420")

# Filter data for enrichment analysis
comparisons.reduced <- lapply(comparisons, function(x) {
  df.filter <- x %>%
    filter(enrichment == "e" & p_fdr_bh < 0.05 & depth > 4) %>%
    select(name, ratio_in_study, ratio_in_pop, depth, study_items, p_fdr_bh)
  df.filter
})

# Combine the filtered data with comparison names
df.list <- mapply(function(x, y) data.frame(cond = rep(x, nrow(y)), y), names, comparisons.reduced, SIMPLIFY = FALSE)
df <- do.call(rbind, df.list)

# DOT plot for GO enrichment
df$name <- factor(df$name, levels = unique(df$name))

bp <- ggplot(df, aes(x = cond, y = name, size = p_fdr_bh, fill = ratio_in_study)) +
  geom_point(shape = 21) +
  scale_size(range = c(5, 1)) +
  labs(x = "", y = "GO TERM") +
  scale_fill_viridis(discrete = FALSE, option = "C") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, family = "serif"),
        axis.text.y = element_text(size = 9, face = "bold", family = "serif"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey80")) +
  scale_y_discrete(limits = rev(levels(df$name)))

# Separate data into White Light and Monochromatic
df.white <- df.list[1:7]
df.w <- do.call(rbind, df.white)
df.mono <- df.list[8:12]
df.m <- do.call(rbind, df.mono)


######### VENN DIAGRAMS - FIG 2C ############
# This section creates Venn diagrams to compare different gene sets

#### 2. Read and filter genes based on log2FoldChange and padj ####

# Read genes for each comparison (470UP, 560UP, 660UP vs 630, 660UP vs 690, PSII control)
genes470_420 <- c470_420$Gene[which(c470_420$log2FoldChange < 0 & c470_420$padj < 0.05)]
genes560_520 <- c560_520$Gene[which(c560_520$log2FoldChange < 0)]
genes660_630 <- c660_630$Gene[which(c660_630$log2FoldChange < 0)]
genes660_690 <- c660_690$Gene[which(c660_690$log2FoldChange < 0)]

# Read PSII genes at different time points (4 min, 30 min, 1 hour)
genes_psii4 <- red_4m$Gene[which(red_4m$log2FoldChange < 0)]
genes_psii30 <- red_30m$Gene[which(red_30m$log2FoldChange < 0)]
genes_psii1h <- red_1h$Gene[which(red_1h$log2FoldChange < 0)]

# Additional gene comparisons without LFC limits
genes_psi30_psi4 <- psi30_psi4$rownames.resSig.[which(abs(psi30_psi4$log2FoldChange) > 1)]
genes_psi60_psi4 <- psi60_psi4$rownames.resSig.[which(abs(psi60_psi4$log2FoldChange) > 1)]
genes_psii30_psii4 <- psii30_psii4$rownames.resSig.[which(abs(psii30_psii4$log2FoldChange) > 1)]
genes_psii60_psii4 <- psii60_psii4$rownames.resSig.[which(abs(psii60_psii4$log2FoldChange) > 1)]

# Compare PSI and PSII 30-min treatments
genes_psi30_psii30 <- psi30_psii30$Gene[which(abs(psi30_psii30$log2FoldChange) > 1)]

# Core gene set and comparison genes
genes420 <- c420_470_660$rownames.resSig.[which(abs(c420_470_660$log2FoldChange) > 1)]
genes_core <- core
genes_contr <- ox_red_c$Gene[which(ox_red_c$log2FoldChange > 0)]

#### 3. Venn Diagram Construction ####

# Create a list of gene sets for comparison
venn_data <- list(genes_psi60_psi4, genes_psii60_psii4, genes420, genes_core)
category_names <- c("psi60_psi4", "psii60_psii4", "c420_470_660", "core")

# Plot the Venn diagram using ggVennDiagram
ggVennDiagram(venn_data, label_alpha = 0, label = "count", category.names = category_names) +
  ggplot2::scale_fill_gradient(low = "white", high = "red")

#### 4. Literature-based Venn Diagram ####

# Load gene lists from literature comparison
piippo <- read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 2)[,1]
fey <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 3)[,1])
bode <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 4)[,1])
adamiec <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 5)[,1])
rossel2002 <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 6)[,1])
rossel2007 <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 7)[,1])
pesaresi <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 8)[,1])
dietzel <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 9)[,1])
brautigam <- toupper(read.xlsx("Comparison of PSIIw_PSIw_many_papers (1).xlsx", sheet = "Vertailu", cols = 10)[,1])

# Venn diagram data and category names
venn_data_literature <- list(score12red, dietzel, piippo, fey, bode)
category_names_literature <- c("psi", "Dietzel", "Piippo", "Fey", "Bode")

# Generate the Venn diagram from literature overlaps
VennDiagram::venn.diagram(x = venn_data_literature, filename = 'round2/MC_modeling/psi_ven.png', category.names = category_names_literature)

#### 5. Enhanced Venn Plot with venneuler ####

# Prepare data for venneuler visualization
Hmat <- lapply(seq_along(venn_data), function(i) {
  Set <- category_names[i]
  data.frame(Element = venn_data[[i]], Set)
})
Hmat <- do.call(rbind, Hmat)

# Generate and plot the venneuler diagram
vd <- venneuler(Hmat)
plot(vd)

# Annotate the labels with set sizes
HH <- Hmat %>%
  group_by(Set) %>%
  summarise(n = n()) %>%
  .[match(vd$labels, .$Set),] %>%
  mutate(n = paste(Set, "\n", n))

# Update plot labels and display the plot
vd$labels <- HH$n
plot(vd)

#### 6. Final Venn Diagram for Manuscript ####

# Final Venn diagram comparing various gene sets for the manuscript figure
VennDiagram::venn.diagram(
  x = list(genes_psi30_psi4, genes_psi60_psi4, genes_psii30_psii4, genes_psii60_psii4, core, genes420),
  filename = 'manuscript/figures/2C_4min.png',
  output = TRUE,
  category.names = c("Set 1", "Set 2", "Set 3", "Set 4", "Set 5", "Set 6")
)

###
##heatmap 
p <- ggplot(moduletrait.df, aes(x=Cluster, y=PQ_time, fill=Cor.))+
  geom_tile(colour="white", size =0.5, width=0.8)+
  guides(fill=guide_legend(title=""))+
  labs(x="", y="", title="")+
  scale_y_discrete(expand=c(0, 0))+
  geom_text(aes(label = round(Cor., 2)), color = "white",size=6) +
  scale_fill_gradient2(low="blue",mid="white", high="red")+
  #coord_fixed()+
  theme_grey(base_size=10)+
  theme(legend.position="none",
        #, legend.direction="vertical",
        #  legend.title=element_text(colour=textcol),
        #  legend.margin=margin(grid::unit(0, "cm")),
        # legend.text=element_text(colour=textcol, size=7, face="bold"),
        #legend.key.height=grid::unit(0.8, "cm"),
        #legend.key.width=grid::unit(0.2, "cm"),
        axis.text.x=element_text(size=10, colour=textcol),
        axis.text.y=element_text(size=10, vjust=0.2, colour=textcol),
        axis.ticks=element_line(linewidth =0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        #plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm" ),
        plot.title=element_text(colour=textcol, hjust=0, size=14, face="bold")
  )

p <- ggplot(moduletrait.df, aes(x = Cluster, y = PQ_time, fill = Cor.)) +
  geom_tile(colour = "white", size = 0.2, width = 0.9, height = 0.2) +  # Adjust tile size
  geom_text(aes(label = round(Cor., 2)), color = "black", size = 3) +  # Adjust text size
  scale_fill_gradient2(low = rgb(0, 0.4, 0.9), mid = "white", high = "red") +
  #scale_fill_viridis(option = "plasma", name = "Correlation", direction = -1) +
  scale_y_discrete(expand = c(0, 0)) +
  guides(fill = guide_legend(title = "")) +
  labs(x = "", y = "", title = "") +
  theme_grey(base_size = 10) +
  theme(
    #plot.margin = unit(c(0,0,0, 0),"inches"), 
    plot.margin = unit(c(0, 1, -1, 1), "cm"),
    axis.title.y = element_text(angle = 0) ,
    legend.position = "none",
    axis.text.x = element_text(size = 12, colour = "black"),  # Rotate x-axis text
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(linewidth = 0.4),
    plot.background = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(colour = "black", hjust = 0, size = 14, face = "bold")
  )



