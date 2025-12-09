#' Fig6ab_rain_box
#' 2025-09-11
#' Qingwang Chen

# setwd
getwd()
setwd("/vast/projects/quartet_rna_refdata/analysis/")

# Specify the new library path
custom_lib_path <- "/vast/projects/quartet_rna_refdata/my_r_packages"

# Modify .libPaths to prioritize the custom library path
.libPaths(c(custom_lib_path, .libPaths()))

# Print the current library paths
print(.libPaths())

library(devtools)
devtools::install_github('erocoar/gghalves')
library(gghalves)

df <- read.table("/vast/projects/quartet_rna_refdata/analysis/R/figures/fig6/data/fig6_application_all_metrics_RefData.csv",
                 sep=",",header = T)
head(df)
fill_group <- c(LO="#4C72B0", LS="#00A087", SO="#DD8452")
proto_levels <- c("D_ONT","M_PAB","P_ONT","I_PAB","P_ILM","P_BGI","R_ILM","R_BGI","R_ELE")
table(df$protocols_platforms,df$group)
proto_shapes <- setNames(c(16,17,15,18,0,1,2,3,4), proto_levels)

df <- df[!df$Metric=="FNR",]
df$Metric <- factor(df$Metric,levels = c("RC","RMSE","MCC"))
df$type <- factor(df$type,levels = c("isoform","AS"))
df$protocols_platforms <- factor(df$protocols_platforms,levels = c("D_ONT","M_PAB","P_ONT","I_PAB","P_ILM","P_BGI","R_ILM","R_BGI","R_ELE"))
pal_simpsons()(16)[c(1:16)]

p <- ggplot(data = df,       
            aes(x=group, y=value, fill=group)) +
  #geom_half_violin(side="r")+
  #geom_half_boxplot(side = "r") +    
  # geom_half_point_panel(side = "l")+
  geom_half_violin(side = "r", color=NA, alpha=0.35) +    
  geom_half_boxplot(side = "r", errorbar.draw = FALSE, 
                    outlier.size = 1,
                    outlier.stroke = 0.2,
                    width=0.2, linewidth=0.5) +    
  geom_half_point_panel(aes( fill = protocols_platforms),side = "l", 
                        range_scale = .85,
                        shape=21, size=1.5, color="white")+
  scale_fill_manual(values = c(pal_simpsons()(16)[c(1,3,2,5,10,7,4,11,16)],#,13,4,15,16
                               c("#4C72B0","#00A087", "#DD8452"))) +
  theme_bw()+
  theme(
    axis.text.x = element_text(face = "bold",color = "black"),
    axis.text.y = element_text(face = "bold",color = "black"),
    axis.title  = element_text(face = "bold",color = "black"),
    panel.grid.minor  =  element_blank(),
    legend.background =  element_blank(),
    legend.title = element_text(face = "bold",color = "black"),
    legend.text = element_text(color = "black")
  )+
  facet_wrap(type~Metric,scales = "free");p
ggsave("/Users/xr/Desktop/Quartet_as/new/Fig6_rain_box_0711.png",p,width = 9,height = 5.5)
topptx(p,"/Users/xr/Desktop/Quartet_as/new/Fig6_rain_box_0711.pptx",width = 7.2,height = 4.4)

