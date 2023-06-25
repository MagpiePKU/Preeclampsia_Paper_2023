
library(Seurat) 
library(ArchR) 
library(dplyr) 
library(ggplot2) 
library(ggthemes) 
library(reshape2) 
library(ggsignif)
library(RColorBrewer)
library(gghalves)
setwd("./Figure2/") 

# Figure 2a 
tropho_metadata <- readRDS("./tropho_metadata.rds")
nb.cols=length(unique(tropho_metadata$seurat_clusters))
mycolors_groups <- colorRampPalette(brewer.pal(name = "Set1",n=9))(nb.cols)
names(mycolors_groups) <-  tropho_metadata$seurat_clusters %>% unique 

pdf('./Figure2a.pdf',width=7,height=6)
ggplot(tropho_metadata, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 0.25, aes(color = seurat_clusters)) + theme_classic() + scale_color_manual(values = mycolors_groups)
dev.off()

# Figure 2b 
tropho_distribution_df <- readRDS("./CTRL_Tropho_latentime_distribution_model.rds")

pdf("./Figure2b.pdf", width = 6, height = 5)
ggplot(tropho_metadata) + 
  geom_histogram(aes(x=latent_time,fill=new_PE,group=new_PE,y=..density..),position='dodge',color='black',size=0.5,alpha=0.6) +
  geom_line(data= tropho_distribution_df,mapping = aes(x=x,y=y, color = distribution, group = distribution),size=1.5) + 
  scale_fill_manual(values = c("CTRL" = "black", "PE" = "red")) + scale_color_manual(values = c('Immature' = "cornflowerblue", 'Mature' = "orange")) + theme_classic() + theme(text = element_text(size = 20),legend.position = "bottom", legend.direction = "vertical") + ylab('Distribution Density') + xlab('Latent Time')  
dev.off()

# Figure 2c 
factor(tropho_metadata$new_cluster, levels = c('VCTp_4','VCTp_8','VCTp_1','VCTp_14','VCT_2', 'VCT_10','VCT_5','EVT_11','VCT_7','VCT_0','VCT_3','SCT_6','EVT_13','EVT_9','SCT_12'), ordered = T) -> tropho_metadata$new_cluster
my_color_set = c("VCT_2" = "#E41A1C",
                 "VCTp_8" = "#815375",
                 "VCT_0" = "#3A85A8",
                 "VCTp_1" = "#46A169",
                 "VCT_3" = "#629363",
                 "VCT_7" = "#8D5B96",
                 "EVT_13" = "#C4625D",
                 "SCT_6" = "#FF7F00",
                 "VCTp_4" = "#FFC81D",
                 "EVT_9" = "#F2E631",
                 "EVT_11" = "#BF862B",
                 "VCT_10" = "#BD6253",
                 "VCT_5" = "#EB7AA9",
                 "VCTp_14" = "#CE8BAE",
                 "SCT_12" = "#999999"
)


pdf("./Figure2c.pdf", height = 5, width = 8)
ggplot() + 
  geom_half_violin(data = tropho_metadata[tropho_metadata$new_PE %in% "CTRL", ], side = "l", scale = "width", aes(fill = new_cluster,x = new_cluster, y = CytoTRACE, alpha = new_PE, color = new_cluster)) + 
  geom_half_violin(data = tropho_metadata[tropho_metadata$new_PE %in% "PE", ], side = "r", scale = "width", aes(fill = new_cluster,x = new_cluster, y = CytoTRACE, alpha = new_PE, color = new_cluster)) + 
  theme_bw() + scale_fill_manual(values = c(my_color_set)) + 
  scale_color_manual(values = my_color_set) + 
  scale_alpha_manual(values = c("CTRL" = 1, "PE" = 0.4)) + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1), legend.title = element_blank(), plot.margin = unit(c(0.6,0.6,0.6,0.6), "cm"), legend.position = "none", axis.title.x = element_blank())  
dev.off()


# Figure 2d 
readRDS("./Trophoblast_EpiTrace_result.rds") -> tropho_epitrace_res
pdf("./Figure2d.pdf", height = 4, width = 8)
ggplot(tropho_epitrace_res[!is.na(tropho_epitrace_res$EpiTraceAge_iterative),],aes(x = paste(Phenotype), y = EpiTraceAge_iterative/gw_plot)) + geom_violin(scale = "width", size = 0.9, aes(fill = Phenotype), alpha = 0.8) + geom_boxplot(outlier.alpha = 0, width = 0.2, color = "black", fill = "grey", size = 0.8) + facet_wrap(~Celltype, ncol = 4, scales = "free_x") + scale_fill_manual(values = c("dimgrey","tomato")) + theme_classic() + theme(legend.position = "bottom", axis.title.x = element_blank(), axis.text.x = element_blank(), text = element_text(size = 20)) + ylab("Mitosis aging rate") + geom_signif(comparisons = list(c("CTRL","PE")), test = "wilcox.test",map_signif_level = T, size = 1, textsize = 10, tip_length = 0, vjust = 0.6)
dev.off()


# Figure 2e 
ggplot(tropho_metadata[tropho_metadata$source%in%c('BYS','HX') & tropho_metadata$PE == 0,], 
       aes(x = latent_time, y = CytoTRACE)) + 
  scale_color_gradient2(low = "steelblue", mid = "orange", high = 'firebrick4', midpoint = 0.5) + 
  geom_point(alpha=0.4,size=1,aes(color = latent_time)) + 
  theme_classic() + 
  geom_density2d(aes(group=1),color='black',size=0.4, linetype = "dashed") +
  theme(text = element_text(size = 20), legend.position = "none",plot.margin = unit(c(0,0.2,0,0.5),"cm")) + 
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  scale_y_continuous(breaks = seq(0,1,0.1)) -> p7

ggplot(tropho_metadata[tropho_metadata$source%in%c('BYS','HX')& tropho_metadata$PE == 0,], aes(x = latent_time, fill=as.character(PE))) + 
  geom_density(size= 1, fill = "white", aes(color = as.character(PE))) + 
  scale_color_manual(values = c("0" = "black", '1' = "red")) + 
  geom_histogram(aes(y=..density..),position='dodge',color='black',size=0.3,alpha=0.4, binwidth = 0.02) + 
  scale_fill_manual(values=c('0'='dimgrey','1'='tomato2')) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  theme(text = element_text(size = 20), legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),plot.margin = unit(c(0.5,0,0,0),"cm")) -> p8

ggplot(tropho_metadata[tropho_metadata$source%in%c('BYS','HX')& tropho_metadata$PE == 0,], aes(x = CytoTRACE, fill=as.character(PE))) + 
  geom_density(size= 1, fill = "white", aes(color = as.character(PE))) + 
  scale_color_manual(values = c("0" = "black", '1' = "red")) + 
  geom_histogram(aes(y=..density..),position='dodge',color='black',size=0.3,alpha=0.4, binwidth = 0.02) + 
  scale_fill_manual(values=c('0'='dimgrey','1'='tomato2')) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  coord_flip() + 
  theme(text = element_text(size = 20), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(),plot.margin = unit(c(0,0.5,0,0.1),"cm")) -> p9

p5 <- ggplot() + theme_classic() + theme(axis.line = element_blank())

patchwork::wrap_plots(list(p7,p9), ncol = 2, widths = c(2,0.6)) -> p10
patchwork::wrap_plots(list(p8,p5), ncol = 2, widths = c(2,0.6)) -> p11

pdf("./Figure2e.pdf", width = 9, height = 5.5)
patchwork::wrap_plots(list(p11,p10), ncol = 1) + patchwork::plot_layout(heights = c(1,2)) + patchwork::plot_annotation("CTRL")
dev.off() 


# Figure 2f 
ggplot(tropho_metadata[tropho_metadata$source%in%c('BYS','HX') & tropho_metadata$PE == 1,], 
       aes(x = latent_time, y = CytoTRACE)) + 
  scale_color_gradient2(low = "steelblue", mid = "orange", high = 'firebrick4', midpoint = 0.5) + 
  geom_point(alpha=0.4,size=1,aes(color = latent_time)) + 
  theme_classic() + 
  geom_density2d(aes(group=1),color='black',size=0.4, linetype = "dashed") +
  theme(text = element_text(size = 20), legend.position = "none",plot.margin = unit(c(0,0.1,0,0.5),"cm")) + 
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  scale_y_continuous(breaks = seq(0,1,0.1)) -> p1


ggplot(tropho_metadata[tropho_metadata$source%in%c('BYS','HX')& tropho_metadata$PE == 1,], aes(x = latent_time, fill=as.character(PE))) + 
  geom_density(size= 1, fill = "white", aes(color = as.character(PE))) + 
  scale_color_manual(values = c("0" = "black", '1' = "red")) + 
  geom_histogram(aes(y=..density..),position='dodge',color='black',size=0.3,alpha=0.4, binwidth = 0.02) + 
  scale_fill_manual(values=c('0'='dimgrey','1'='tomato2')) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  theme(text = element_text(size = 20), legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(),plot.margin = unit(c(0.5,0,0,0),"cm")) -> p2

ggplot(tropho_metadata[tropho_metadata$source%in%c('BYS','HX')& tropho_metadata$PE == 1,], aes(x = CytoTRACE, fill=as.character(PE))) + 
  geom_density(size= 1, fill = "white", aes(color = as.character(PE))) + 
  scale_color_manual(values = c("0" = "black", '1' = "red")) + 
  geom_histogram(aes(y=..density..),position='dodge',color='black',size=0.3,alpha=0.4, binwidth = 0.02) + 
  scale_fill_manual(values=c('0'='dimgrey','1'='tomato2')) + 
  theme_classic() + 
  scale_x_continuous(breaks = seq(0,1,0.1)) + 
  coord_flip() + 
  theme(text = element_text(size = 20), legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(),plot.margin = unit(c(0,0.5,0,0.1),"cm")) -> p3

patchwork::wrap_plots(list(p1,p3), ncol = 2, widths = c(2,0.6)) -> p4 
p5 <- ggplot() + theme_classic() + theme(axis.line = element_blank())

patchwork::wrap_plots(list(p2,p5), ncol = 2, widths = c(2,0.6)) -> p6 

pdf("./Figure2f.pdf", width = 9, height = 5.5)
patchwork::wrap_plots(list(p6,p4), ncol = 1) + patchwork::plot_layout(heights = c(1,2)) + patchwork::plot_annotation("PE")
dev.off() 


# Figure 2g 
tropho_metadata[tropho_metadata$source %in% c('HX','BYS'),] %>% dplyr::group_by(orig.ident,new_PE) %>% dplyr::summarise(immature=sum(latent_time<0.5 & CytoTRACE>0.2),all=sum(CytoTRACE<2)) -> df_cytotrace_freq
df_cytotrace_freq$immature_freq <- df_cytotrace_freq$immature/df_cytotrace_freq$all
factor(df_cytotrace_freq$new_PE, levels = c("PE","CTRL")) -> df_cytotrace_freq$new_PE 

pdf("./Figure2g.pdf", width = 6, height = 3)
ggplot(df_cytotrace_freq,aes(x=immature_freq,y=new_PE,fill=new_PE)) + 
  geom_boxplot() + 
  ggbeeswarm::geom_beeswarm() + 
  theme_classic() + 
  scale_fill_manual(values=c('CTRL' = 'gray','PE' = 'tomato')) + 
  theme(axis.title.x=element_text(size=15),text=element_text(size=25), axis.title.y = element_blank(), legend.position = "none") + 
  scale_x_continuous(position = "top") + 
  xlab('3rd Trimester Immature Trophoblast Prevalence') + 
  ggpubr::stat_compare_means(comparisons = list(c('CTRL','PE')),size=5,label = 'p.value',method = 't.test')
dev.off() 



