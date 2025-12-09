#' fig4-statistics-of-reference-isoforms
#' Qingwang Chen
#' 2025-06-13
#' 

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())

# install.packages("Rcpp", type = "source", lib=custom_lib_path)
# library(Rcpp)
# library import
library(edgeR)
library(arrow)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(grid) 

# --------------------------- 2. 参数与文件 --------------------------------
files <- tribble(
  ~path,                                                                                                  ~type,   ~dataset, ~skip_filter,
  "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/ref_expr_b7_p4_s84_u_20250516.csv",  "refFC",  "LO", TRUE ,
  "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/lo-e-g/RefData_DEIs_all_isoforms_classified_u_20250522.csv", "refDEI", "LO", FALSE,
  "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/ref_expr_b13_p4_s156_u_20250407.csv",  "refFC",  "SO", TRUE ,
  "/vast/projects/quartet_rna_refdata/analysis/isoforms/ref_data_construction/so-e-g/RefData_DEIs_all_isoforms_classified_u_20250407.csv", "refDEI", "SO", FALSE
)

read_and_count <- function(path, type, dataset, skip_filter = FALSE) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  names(df) <- tolower(names(df))          # ← 代替 janitor::clean_names()
  
  if (!"compare" %in% names(df))
    stop("Missing 'compare' column in: ", path)
  
  if (!skip_filter && type %in% c("refDEI", "refDAS")) {
    if (!"final" %in% names(df))
      stop("Missing 'final' column in: ", path)
    df <- df |>
      dplyr::filter(final != ifelse(type == "refDEI", "non-DEI", "non-DAS"))
  }
  
  df |>
    dplyr::count(compare, name = "count") |>
    dplyr::mutate(type = type, dataset = dataset)
}


raw_df <- purrr::pmap_dfr(files, read_and_count)

# --------------------------- 3. 数据整形 ----------------------------------
plot_df <- raw_df |>
  mutate(
    compare = factor(compare, levels = c("D5/D6","F7/D6","M8/D6")),
    dataset = factor(dataset, levels = c("LO","SO")),
    type    = factor(type,    levels = c("refFC","refDEI"))
  )

# 双轴缩放因子
sf <- with(plot_df,
           max(count[type=="refFC"], na.rm = TRUE) /
             max(count[type=="refDEI"],na.rm = TRUE))

plot_df <- plot_df |> 
  mutate(count_scaled = if_else(type=="refDEI", count*sf, count))

# --------------------------- 4. 绘图 --------------------------------------
nbt_cols <- c("D5/D6"="#1b9e77","F7/D6"="#d95f02","M8/D6"="#7570b3")

p <- ggplot(
  plot_df[plot_df$dataset == "LO", ],
  aes(dataset, count_scaled, fill = compare)
) +
  geom_col(position = position_dodge(.7), width = .55) +
  facet_wrap(
    ~ type, scales = "free_x",
    labeller = labeller(
      type = c(refFC  = "refFC Isoform",
               refDEI = "refDEI")
    )
  ) +
  scale_fill_manual(values = nbt_cols) +
  scale_y_continuous(
    name     = "refFC Isoform Count",
    sec.axis = sec_axis(~ . / sf, name = "refDEI Count")
  ) +
  labs(x = "Dataset", fill = "Sample pair") +   # ← 显示 x 轴标题
  theme_minimal(base_size = 14) +
  theme(
    ## 文字
    axis.text.x  = element_text(size = 12, face = "bold"),
    axis.text.y  = element_text(size = 12),
    strip.text   = element_text(face = "bold"),
    legend.position = "top",
    
    ## 轴线 & 须线
    axis.line.x  = element_line(colour = "black", linewidth = 0.6),
    axis.line.y  = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )

p+canvas(width = 6,height = 4)

ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4b_isoform_num_LO.png",p, width=6, height=4, dpi=300)
ggsave("/vast/projects/quartet_rna_refdata/analysis/figures/fig4/fig4b_isoform_num_LO.pdf",p, width=6, height=4, dpi=300)

# remove all
rm(list=ls())
gc()