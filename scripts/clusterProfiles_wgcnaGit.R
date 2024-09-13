# 
library(ggplot2)
library(gridExtra)
library(viridisLite)
library(reshape2)

# Load data, results from clustering, with genes asssigned to clusters.
load("MEs.RData")

# Reorder and rename MEs row names
order <- c(26:33,38:44, 46:49,34, 45, 50:51, 52:55, 35:37, 56, 3:5, 1,6,7, 8:10, 2, 11:12, 13:25)
MEs <- MEs[order,]
rownames(MEs)[9:31] <- c("PSI_4a", "PSI_4b","PSI_4c", "PSI_4e", "PSI_30a","PSI_30b", "PSI_30c", "PSI_30d",
                         "PSI_60a","PSI_60b", "PSI_60c", "PSII_4a", "PSII_4b","PSII_4c", "PSII_4d", "PSII_30a",
                         "PSII_30b", "PSII_30c", "PSII_30d", "1PSII_60b", "PSII_60c", "PSII_60d", "PSII_60a")

rownames(MEs)[32:56] <- c("420nm-2", "420nm_3","420nm-4","470nm_2", "470nm-1", "470nm-3", "520nm-1", "520nm-",
                          "520nm-3u", "560nm_1", "560nm-2", "560nm-3", "630nm-A1", "630nm-B1", "630nm-D1", 
                          "660nm-1", "660nm-3", "660nm-4", "690nm-1", "690nm-2", "690nm-4", "white-1", 
                          "white-2184", "white-2", "white-4")

groups <- c(rep("OX",4), rep("RED",4), rep("PSI_4",4), rep("PSI_30",4), rep("PSI_60",3), rep("PSII_4",4),
            rep("PSII_30",4), rep("PSII_60",4), rep("420L",3), rep("470L",3), rep("520L",3), rep("560L",3),
            rep("630L",3), rep("660L",3), rep("690L",3), rep("White",4))

# Generate plots
pal <- viridis(36)
g_plots <- list()

for (i in 1:ncol(MEs)) {
  module.df <- data.frame(Samples = rownames(MEs), eigengene = MEs[, i])
  
  # Create spline
  #spline_int <- as.data.frame(spline(1:nrow(module.df), module.df$eigengene))
  
  # Plotting
  bp <- ggplot(module.df) + 
    geom_point(aes(x = Samples, y = eigengene), size = 1) +
    geom_line(data = spline_int, aes(x = x, y = y), color = pal[i], alpha = 0.4) + 
    theme(axis.text.x = element_text(angle = 45, size = 6)) + 
    ggtitle(colnames(MEs)[i])
  
  g_plots[[i]] <- ggplotGrob(bp)
}


# Save cluster profile plots to PDF
pdf("clusterProfiles.pdf")
marrangeGrob(grobs = g_plots, nrow = 3, ncol = 3)
dev.off()
