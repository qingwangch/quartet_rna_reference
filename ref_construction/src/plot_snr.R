#' plot_snr.R
#' 


# mapping relationship
batch_shape_mapping <- c(
  # 已有映射
  "D_ONT_LG_B1" = 16,
  "D_ONT_LW_B1" = 15,
  "I_PAB_LN_B1" = 17,
  "I_PAB_LN_B2" = 18,
  "M_PAB_LG_B1" = 3,
  "M_PAB_LG_B2" = 4,
  "M_PAB_LN_B1" = 8,
  "P_ONT_LB_B1" = 25,
  "P_ONT_LG_B1" = 22,
  "P_ONT_LG_B2" = 23,
  "P_ONT_LN_B1" = 0,
  "P_ONT_LN_B2" = 1,
  "P_BGI_L3_B1" = 2,
  "P_ILM_L8_B1" = 5,
  "R_BGI_L3_B1" = 6,
  "R_ILM_L8_B1" = 7,
  # 新增映射
  "P_BGI_L6_B1" = 9,
  "P_ILM_L1_B1" = 10,
  "P_ILM_L2_B1" = 11,
  "P_ILM_L5_B1" = 12,
  "P_ILM_L6_B1" = 13,
  "R_BGI_L6_B1" = 14,
  "R_BGI_L7_B1" = 19,
  "R_ELE_LH_B1" = 20,
  "R_ILM_L1_B1" = 21,
  "R_ILM_L2_B1" = 24,
  "R_ILM_L2_B2" = 16,   # =25 已用；改用 26 → 26 不在基础 set
  "R_ILM_L4_B1" = 15,  # 26; fallback: 26
  "R_ILM_L4_B2" = 17,  # 27
  "R_ILM_L4_B3" = 18,  # 28
  "R_ILM_L5_B1" = 3,  # 29
  "R_ILM_L5_B2" = 4,  # 30
  "R_ILM_L5_B3" = 8,  # 31
  "R_ILM_L6_B1" = 25,   # 32
  "C_BGI_LS_B1" = 16,
  "C_BGI_LS_B2" = 15
)

plot_snr_quartet <- function(SNR_x,metadata,exprMat,output_prefix,width,height) {
  
  SNR_n <- round(SNR_x$signoise_db, 1)
  pca.all <- SNR_x$pca_prcomp
  pcs <- pca.all$x %>% as.data.frame()
  pcs <- pcs[, 1:3]
  pcs$code <- rownames(pcs)
  pcs <- merge(pcs,metadata,by="code")
  pcs$sample <- factor(pcs$sample, levels=c("D5", "D6", "F7", "M8"))
  
  subtype_pal <- c('#4CC3D9' ,'#7BC8A4' ,'#FFC65D', '#F16745')
  shape_values <- batch_shape_mapping
  
  p <- ggplot(pcs, aes(x = PC1, y = PC2, color = sample)) +
    geom_point(size = 3.5, aes(shape = Batch)) +
    theme_bw() +
    theme(
      axis.text = element_text(),
      legend.position = "right",
      legend.background = element_blank(),
      plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    guides(
      fill = guide_legend("Type", override.aes = list(shape = shape_values, fill = subtype_pal))
    ) +
    scale_color_manual("Type", values = subtype_pal) +
    scale_shape_manual("Batch", values = shape_values) +
    labs(
      x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2, 1] * 100),
      y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2, 2] * 100),
      title = paste0(output_prefix),
      subtitle = paste0("N = ", nrow(exprMat), "  SNR = ", SNR_n)
    ) +
    theme(aspect.ratio = 1 / 1)
  
  p_n <- p +
    theme(
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12)
    )
  
  ggsave(paste0("./results/plots/", output_prefix, ".png"), p, width = width, height = height, dpi = 300)
  
  return(list(p = p_n, pcs = pcs))
}

plot_snr_all <- function(SNR_x,metadata,exprMat,output_prefix,width,height) {
  
  SNR_n <- round(SNR_x$signoise_db, 1)
  pca.all <- SNR_x$pca_prcomp
  pcs <- pca.all$x %>% as.data.frame()
  pcs <- pcs[, 1:3]
  pcs$code <- rownames(pcs)
  pcs <- merge(pcs,metadata,by="code")
  pcs$group <- factor(pcs$group,levels = c("Quartet","HCC1395","HCC1395BL"))
  # color
  subtype_pal <- c("#4091c0","#69001f","#f7b293")
  shape_values <- batch_shape_mapping
  
  p <- ggplot(pcs, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3.5, aes(shape = Batch)) +
    theme_bw() +
    theme(
      axis.text = element_text(),
      legend.position = "right",
      legend.background = element_blank(),
      plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    guides(
      fill = guide_legend("group", override.aes = list(shape = shape_values, fill = subtype_pal))
    ) +
    scale_color_manual("group", values = subtype_pal) +
    scale_shape_manual("Batch", values = shape_values) +
    labs(
      x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2, 1] * 100),
      y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2, 2] * 100),
      title = paste0(output_prefix),
      subtitle = paste0("N = ", nrow(exprMat), "  SNR = ", SNR_n)
    ) +
    theme(aspect.ratio = 1 / 1)
  
  p_n <- p +
    theme(
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12)
    )
  
  ggsave(paste0("./results/plots/", output_prefix, ".png"), p, width = width, height = height, dpi = 300)
  
  return(list(p = p_n, pcs = pcs))
}

plot_snr_hcc1395 <- function(SNR_x,metadata,exprMat,output_prefix,width,height) {
  
  SNR_n <- round(SNR_x$signoise_db, 1)
  pca.all <- SNR_x$pca_prcomp
  pcs <- pca.all$x %>% as.data.frame()
  pcs <- pcs[, 1:3]
  pcs$code <- rownames(pcs)
  pcs <- merge(pcs,metadata,by="code")
  pcs$group <- factor(pcs$group,levels = c("Quartet","HCC1395","HCC1395BL"))
  # color
  subtype_pal <- c("#69001f","#f7b293")
  shape_values <- batch_shape_mapping
  
  p <- ggplot(pcs, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3.5, aes(shape = Batch)) +
    theme_bw() +
    theme(
      axis.text = element_text(),
      legend.position = "right",
      legend.background = element_blank(),
      plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    guides(
      fill = guide_legend("group", override.aes = list(shape = shape_values, fill = subtype_pal))
    ) +
    scale_color_manual("group", values = subtype_pal) +
    scale_shape_manual("Batch", values = shape_values) +
    labs(
      x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2, 1] * 100),
      y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2, 2] * 100),
      title = paste0(output_prefix),
      subtitle = paste0("N = ", nrow(exprMat), "  SNR = ", SNR_n)
    ) +
    theme(aspect.ratio = 1 / 1)
  
  p_n <- p +
    theme(
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12)
    )
  
  ggsave(paste0("./results/plots/", output_prefix, ".png"), p, width = width, height = height, dpi = 300)
  
  return(list(p = p_n, pcs = pcs))
}

plot_snr_six <- function(SNR_x,metadata,exprMat,output_prefix,width,height) {
  
  SNR_n <- round(SNR_x$signoise_db, 1)
  pca.all <- SNR_x$pca_prcomp
  pcs <- pca.all$x %>% as.data.frame()
  pcs <- pcs[, 1:3]
  pcs$code <- rownames(pcs)
  pcs <- merge(pcs,metadata,by="code")
  pcs$sample <- factor(pcs$sample,levels = c("D5", "D6", "F7", "M8", "HCC1395","HCC1395BL"))
  # color
  subtype_pal <- c('#4CC3D9' ,'#7BC8A4' ,'#FFC65D', '#F16745', "#69001f", "#f7b293")
  shape_values <- batch_shape_mapping
  
  p <- ggplot(pcs, aes(x = PC1, y = PC2, color = sample)) +
    geom_point(size = 3.5, aes(shape = Batch)) +
    theme_bw() +
    theme(
      axis.text = element_text(),
      legend.position = "right",
      legend.background = element_blank(),
      plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    guides(
      fill = guide_legend("group", override.aes = list(shape = shape_values, fill = subtype_pal))
    ) +
    scale_color_manual("group", values = subtype_pal) +
    scale_shape_manual("Batch", values = shape_values) +
    labs(
      x = sprintf("PC1 (%.2f%%)", summary(pca.all)$importance[2, 1] * 100),
      y = sprintf("PC2 (%.2f%%)", summary(pca.all)$importance[2, 2] * 100),
      title = paste0(output_prefix),
      subtitle = paste0("N = ", nrow(exprMat), "  SNR = ", SNR_n)
    ) +
    theme(aspect.ratio = 1 / 1)
  
  p_n <- p +
    theme(
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12)
    )
  
  ggsave(paste0("./results/plots/", output_prefix, ".png"), p, width = width, height = height, dpi = 300)
  
  return(list(p = p_n, pcs = pcs))
}
