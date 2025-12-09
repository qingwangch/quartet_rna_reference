calSNR_pca <- function(exprMat, group) {
  library(data.table)
  
  IDs <- colnames(exprMat)
  IDs.group.mat <- data.table(IDs = IDs, group = group)
  
  # Remove features which variance is zero to ensure the PCA
  exprMat <- exprMat[which(apply(exprMat, 1, var) != 0), ]
  
  # PCA
  pca_prcomp <- prcomp(t(exprMat), retx = TRUE, center = TRUE, scale. = FALSE)
  # pca_prcomp <- prcomp(t(exprMat),scale. = T)
  pcs <- as.data.frame(predict(pca_prcomp))
  pcs$Sample_id <- rownames(pcs)
  
  # Percent: Proportion of Variance, AccumPercent: Cumulative Proportion
  dt.perc.pcs <- data.table(PCX = 1:length(summary(pca_prcomp)$importance[2,]),
                            Percent = summary(pca_prcomp)$importance[2,],
                            AccumPercent = summary(pca_prcomp)$importance[3,])
  
  dt.dist <- data.table(ID.A = rep(IDs, each = length(IDs)),
                        ID.B = rep(IDs, time = length(IDs)))
  
  dt.dist$group.A <- IDs.group.mat[match(dt.dist$ID.A, IDs.group.mat$IDs)]$group
  dt.dist$group.B <- IDs.group.mat[match(dt.dist$ID.B, IDs.group.mat$IDs)]$group
  
  dt.dist[, Type := ifelse(ID.A == ID.B, 'Same',
                           ifelse(group.A == group.B, 'Intra', 'Inter'))]
  
  dt.dist[, Dist := (dt.perc.pcs[1]$Percent * (pcs[ID.A, 1] - pcs[ID.B, 1])^2 +
                       dt.perc.pcs[2]$Percent * (pcs[ID.A, 2] - pcs[ID.B, 2])^2)]
  
  dt.dist.stats <- dt.dist[, .(Avg.Dist = mean(Dist)), by = .(Type)]
  setkey(dt.dist.stats, Type)
  signoise <- dt.dist.stats['Inter']$Avg.Dist / dt.dist.stats['Intra']$Avg.Dist
  
  signoise_db <- 10 * log10(signoise)
  
  return(list(signoise_db = signoise_db, pca_prcomp = pca_prcomp))
}
