##Assessing the impact of the experimental setup and other variables on the gene expression data.
##Linear model of the gene expression as a linear combination of the measuremetns: 
##reduction of the the Plastoquinone pool, meased at 4min after illumination, PQ_4m
##Reduction of the the Plastoquinone pool, meased at 1h after illumination, PQ_1h
##Ration of the two forms of phytochrome, Pr and Pfr, Pfr_div_Pr. 
##Concentartion of the blue right receptor Phototropin and Cry2. 
##Electron transfer rate, ERT.

library(plyr)
library(readr)
library(caret)
library(ggplot2)
library(repr)
library(dplyr)
library(glmnet)
library(fitdistrplus)
library(ridge)
library(gridExtra)
library(jtools) 
library(MASS)

##read the TPM normalized expression counts. Expression is averaged over the replicates. 
load("tpm_gene.RData")
expr.avg=tpm.gene
rownames(expr.avg)=expr.avg$gene
expr.avg=expr.avg[,-1]
ids =which(rowSums(expr.avg)==0)
expr.avg=expr.avg[-ids,]

#groups=c("470", "560", rep("420",3), rep("470",2), rep("520",3), rep("560",2), rep("630",3), rep("660",3), rep("690",3), rep("white",4),
        # rep("OX",4),rep("RED",4), "PSII4", rep("PSII60",3), rep("PSI4",4), rep("PSI30",3), "PSII4", "PSI30", rep("PSI60",3), rep("PSII4",2), rep("PSII30",4), "PSII60")

#grp.ord = order(groups)
#expr = data.frame(cbind(grp=groups[grp.ord], datExpr[grp.ord,]))


##load the experimental data from a table.
data = openxlsx::read.xlsx("ForLinearModeling.xlsx", rows = 2:29)
colnames(data)[2]="Pfr_div_Pr"
data=data[1:7,]

# Data Exploration
for (j in ids_pq4) {
  
  cluster.model$cluster1_Expr_A <- t(expr.avg[j, 1:9])
  
  # Outliers
  par(mfrow = c(1, 2))
  boxplot(cluster.model$cluster1_Expr_A, 
          xlab = "", 
          ylab = "Expression",
          main = "Boxplot of Expression")
  dotchart(cluster.model$cluster1_Expr_A,
           xlab = "Expression",
           ylab = "Order of Data",
           main = "Dotchart of Expression")
  
  # Relationship
  p <- ggplot(data = cluster.model, aes(x = PQ_4m, y = cluster1_Expr_A)) +
    geom_point() +
    xlab("PQ_4m") + ylab("Expression") +
    ggtitle("Relationship between PQ_4m and Expression") +
    theme(text = element_text(size = 15))
  print(p)
  
  # Interactions
  p <- ggplot(data = data, aes(x = PQ_1h, y = Pfr_div_Pr)) +
    geom_point() +
    xlab("PQ_1h") + ylab("Pfr_div_Pr") +
    ggtitle("Interaction between PQ_1h and Pfr_div_Pr") +
    geom_smooth(method = "lm") +
    theme(text = element_text(size = 15))
  print(p)
}

# Model Validation
for (i in 1:length(ids_1h)) {
  print(i)
  cluster.model$cluster1_Expr_A <- t(expr.avg[ids_1h[i], 1:9])
  
  # Fit model
  M1 <- lm(cluster1_Expr_A ~ PQ_4m, data = cluster.model)
  E1 <- rstandard(M1)
  F1 <- fitted(M1)
  
  # Homogeneity
  plot(F1, E1, 
       xlab = "Fitted Values", 
       ylab = "Residuals", 
       main = "Residuals vs Fitted Values")
  
  # Normality
  hist(E1, 
       main = "Histogram of Residuals", 
       xlab = "Residuals", 
       breaks = 10)
  
  # Dependency
  plot(cluster.model$PQ_4m, E1, 
       xlab = "PQ_4m", 
       ylab = "Residuals", 
       main = "Residuals vs PQ_4m")
  abline(h = 0, lty = 2)
}
####

# Adjusted R square data frames
adj.rsquare <- data.frame(
  model_Pfr_div_Pr = numeric(),
  mod_Cry2 = numeric(),
  mod_PQ_4m = numeric(),
  mod_PQ_1h = numeric(),
  mod_etr = numeric()
)

adj.rsquare.comb <- data.frame(
  pq_etr = numeric(),
  pq_cry = numeric(),
  pq_pfr = numeric()
)

# Iterate through each row of gene expression data, expr.avg
for (i in 1:nrow(expr.avg)) {
  
  print(i)
  
  # Prepare expression data for current row
  cluster.model$cluster1_Expr_A <- t(expr.avg[i, 1:9])
  
  # Fit models with single predictors
  mod_Phototropin <- lm(cluster1_Expr_A ~ Phototropin, data = cluster.model)
  mod_LeafAbs <- lm(cluster1_Expr_A ~ LeafAbs, data = cluster.model)
  mod_Pfr_div_Pr <- lm(cluster1_Expr_A ~ Pfr_div_Pr, data = cluster.model)
  mod_Cry2 <- lm(cluster1_Expr_A ~ Cry2, data = cluster.model)
  mod_PQ_4m <- lm(cluster1_Expr_A ~ PQ_4m, data = cluster.model)
  mod_PQ_1h <- lm(cluster1_Expr_A ~ PQ_1h, data = cluster.model)
  mod_etr <- lm(cluster1_Expr_A ~ ETR, data = cluster.model)
  
  # Extract adjusted R-squared values and p-values
  summary_PQ_1h <- summary(mod_PQ_1h)
  pq1h_pval <- summary_PQ_1h$coefficients[2, 4]
  r2_PQ_1h <- summary_PQ_1h$r.squared
  
  # Check if the model with PQ_1h is significant and explains enough variation
  if (!is.na(r2_PQ_1h) && r2_PQ_1h > 0.7 && pq1h_pval < 0.05) {
    
    # Fit models with PQ_1h and additional predictors
    mod_pq1h_etr <- lm(cluster1_Expr_A ~ ETR + PQ_1h, data = cluster.model)
    mod_pq1h_cry <- lm(cluster1_Expr_A ~ Cry2 + PQ_1h, data = cluster.model)
    mod_pq1h_pfr <- lm(cluster1_Expr_A ~ Pfr_div_Pr + PQ_1h, data = cluster.model)
    
    # Extract adjusted R-squared values
    r2_pq1h_etr <- summary(mod_pq1h_etr)$r.squared
    r2_pq1h_cry <- summary(mod_pq1h_cry)$r.squared
    r2_pq1h_pfr <- summary(mod_pq1h_pfr)$r.squared
    
    # Update adj.rsquare.comb if any combination improves upon the single predictor model
    if (max(r2_pq1h_etr, r2_pq1h_cry, r2_pq1h_pfr) > r2_PQ_1h) {
      adj.rsquare.comb[i, ] <- c(r2_pq1h_etr, r2_pq1h_cry, r2_pq1h_pfr)
    } else {
      adj.rsquare.comb[i, ] <- c(0, 0, 0)
    }
    
  } else {
    adj.rsquare.comb[i, ] <- c(0, 0, 0)
  }
  
  # Store adjusted R-squared values for single predictor models
  adj.rsquare[i, ] <- c(
    summary(mod_Pfr_div_Pr)$r.squared,
    summary(mod_Cry2)$r.squared,
    summary(mod_PQ_4m)$r.squared,
    summary(mod_PQ_1h)$r.squared,
    summary(mod_etr)$r.squared
  )
}

# Define threshold
threshold <- 0.7

# Identify rows where models meet the R-squared threshold
ids_pfr <- which(adj.rsquare$model_Pfr_div_Pr > threshold)
ids_cry <- which(adj.rsquare$mod_Cry2 > threshold)
ids_pq4 <- which(adj.rsquare$mod_PQ_4m > threshold)
ids_1h <- which(adj.rsquare$mod_PQ_1h > threshold)
ids_photot <- which(adj.rsquare$photoptr > threshold)
ids_leafAbs <- which(adj.rsquare$leafAbs > threshold)

# Identify the best combination of two covariates
adj.rsquare.max <- apply(adj.rsquare.comb, 1, max, na.rm = TRUE)
max.id <- apply(adj.rsquare.comb, 1, which.max)

etr <- which(max.id == 1 & adj.rsquare.max != 0)
cry <- which(max.id == 2 & adj.rsquare.max != 0)
pfr <- which(max.id == 3 & adj.rsquare.max != 0)

# Compute adjusted R-squared for a specific combination of covariates
adj.rsquare3 <- numeric()
for (id in etr) {
  
  print(id)
  
  cluster.model$cluster1_Expr_A <- t(expr.avg[id, 1:7])
  adj.rsquare3 <- c(adj.rsquare3, summary(lm(cluster1_Expr_A ~ ETR + PQ_1h + LeafAbs, data = cluster.model))$r.squared)
}



##GLM, fitting negative binomial models of the gene expression data.
# Initialize vectors for p-values and minimum p-values
pvalues <- rep(1, nrow(expr.avg))
min.pval <- rep(0, nrow(expr.avg))

# Iterate through each row of expr.avg
for (i in 1:nrow(expr.avg)) {
  
  print(i)
  
  # Prepare expression data for the current row
  cluster.model$cluster1_Expr_A <- round(expr.avg[i, 1:9])
  
  skip_to_next <- FALSE
  
  # Try to fit negative binomial models
  mod_Pfr_div_Pr <- tryCatch(glm.nb(cluster1_Expr_A ~ Pfr_div_Pr, data = cluster.model), 
                             error = function(e) { skip_to_next <<- TRUE })
  mod_Cry2 <- tryCatch(glm.nb(cluster1_Expr_A ~ Cry2, data = cluster.model), 
                       error = function(e) { skip_to_next <<- TRUE })
  mod_PQ_4m <- tryCatch(glm.nb(cluster1_Expr_A ~ PQ_4m, data = cluster.model), 
                        error = function(e) { skip_to_next <<- TRUE })
  mod_PQ_1h <- tryCatch(glm.nb(cluster1_Expr_A ~ PQ_1h, data = cluster.model), 
                        error = function(e) { skip_to_next <<- TRUE })
  
  if (skip_to_next) { next }
  
  # Extract p-values and find the minimum p-value and its index
  pvals <- c(
    summary(mod_Pfr_div_Pr)$coefficients[2, 4],
    summary(mod_Cry2)$coefficients[2, 4],
    summary(mod_PQ_4m)$coefficients[2, 4],
    summary(mod_PQ_1h)$coefficients[2, 4]
  )
  
  pvalues[i] <- min(pvals)
  min.pval[i] <- which.min(pvals)
}

# Identify rows with significant p-values (threshold < 0.05) and the smallest p-value comes from PQ_1h (index 4)
pq1.nb <- rownames(expr.avg)[which(pvalues < 0.05 & min.pval == 4)]
