

library(dplyr) 
library(ggplot2) 
library(ArchR) 
library(patchwork)

setwd("./Figure7/") 


# Figure 7a 
require(reticulate)
np <- import("numpy")
scv <- import('scvelo')

scv$read('/gpfs/output/SingleCell_Placenta/analyse/Merge_Tropho_and_Endo_for_Paper/Save/Arterial_CTRL_and_PE_20220523.h5ad') -> combined_adata
color_list_a = c( "dimgrey", "tomato","red","grey","lightgrey")
output_name <- 'Figure7a'

scv$pl$velocity_embedding_stream(combined_adata, basis='umap',dpi=1000,save=paste0(output_name,'.UMAP.combined.Streamline.reticulate.png'),linewidth=0.5,legend_loc='on data', color = "annotation", color_map='viridis', size = 18, alpha = 1, palette = color_list_a,density = 0.8)

# Figure 7b 
arterial_scv_meta_df <- readRDS("./Arterial_Velocity_metadata.rds")

c( 'Arterial_ACE'="dimgrey", 'Arterial_ACE_PAPPA2'="tomato",'Arterial_ACE_PAPPA2_CLIC3'="red",'Arterial_BMP6'="grey",'Arterial_CD24'="lightgrey") -> MyEndoColorPalette

ggplot(arterial_scv_meta_df[arterial_scv_meta_df$source %in% c("BYS","HX"), ], aes(x = annotation, y = root_cells)) + geom_violin(scale = "width", size = 1, aes(fill = annotation)) + geom_boxplot(width = 0.25, size = 0.5, fill = "black", outlier.alpha = 0, color = "white") + theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), text = element_text(size = 20), legend.position = "none") + scale_fill_manual(values = MyEndoColorPalette) + ylab("Root Prob.") -> r1

pdf("./Figure7b.pdf", width = 6, height = 6)
r1
dev.off()

# Figure 7c 
ggplot(arterial_scv_meta_df[arterial_scv_meta_df$source %in% c("BYS","HX"), ], aes(x = annotation, y = end_points)) + geom_violin(scale = "width", size = 1, aes(fill = annotation)) + geom_boxplot(width = 0.25, size = 0.5, fill = "black", outlier.alpha = 0, color = "white") + theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), text = element_text(size = 20), legend.position = "none") + scale_fill_manual(values = MyEndoColorPalette) + ylab("Terminal Prob.") -> t1

pdf("./Figure7c.pdf", width = 6, height = 6)
t1
dev.off()

# Figure 7d 
table(arterial_scv_meta_df$new_PE, arterial_scv_meta_df$annotation) %>% as.data.frame -> endo_summary_df 
colnames(endo_summary_df) <- c("Phenotype", "annotation", 'Number') 
factor(endo_summary_df$Phenotype, levels = c("CTRL","PE")) -> endo_summary_df$Phenotype 
factor(endo_summary_df$annotation, levels = c("Arterial_CD24",'Arterial_BMP6', "Arterial_ACE","Arterial_ACE_PAPPA2", "Arterial_ACE_PAPPA2_CLIC3")) -> endo_summary_df$annotation

pdf("./Figure7d.pdf", width = 6, height = 6)
ggplot(endo_summary_df, aes(x = annotation, y = Number)) + geom_col(color = "black", aes(fill = Phenotype), size = 0.7, position = "dodge2",width = 0.9, alpha = 0.8)  + theme_classic() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), text = element_text(size = 20), legend.position = "bottom") + scale_fill_manual(values = c("CTRL" = "lightgrey", "PE" = "tomato")) + ylab("Cell Number")
dev.off() 

# Figure 7e 

load('./ACE_pos_Arterial_GSEA_result_20220606.RData')
HallMarkRes_PE_Art_ACE_sig[abs(HallMarkRes_PE_Art_ACE_sig$NES) > 1 & HallMarkRes_PE_Art_ACE_sig$pval < 0.05, ] -> HallMarkRes_PE_Art_ACE_sig_pick
HallMarkRes_PE_Art_ACE_sig_pick$label = gsub("HALLMARK_","", HallMarkRes_PE_Art_ACE_sig_pick$pathway)
HallMarkRes_PE_Art_ACE_sig_pick %>% arrange(NES) -> HallMarkRes_PE_Art_ACE_sig_pick 
HallMarkRes_PE_Art_ACE_sig_pick$label <- factor(HallMarkRes_PE_Art_ACE_sig_pick$label, levels = rev(HallMarkRes_PE_Art_ACE_sig_pick$label))

pdf("./Figure7e.pdf", width = 24, height = 10)
ggplot(HallMarkRes_PE_Art_ACE_sig_pick, aes(x = NES, y = label)) + 
  geom_col(color = "black", size = 0.8, aes(fill = -sign(NES)*log10(pval))) + 
  # geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "lightgrey", size = 1.5) + 
  theme_classic() + scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "red", midpoint = 0) + 
  theme(text = element_text(size = 35), legend.position = "right", legend.direction = "vertical") + 
  ylab("MSigDB Hallmark Pathway") + ggtitle("Arteriole Ec GSEA")
dev.off()


# Figure 7h 
PE_vs_CTRL_art_exp_df <- readRDS("./PE_vs_CTRL_art_exp_df_DE.rds") 
AVM_vs_CTRL_art_exp_df <- readRDS("./AVM_vs_CTRL_art_exp_df_DE.rds")

ggplot(PE_vs_CTRL_art_exp_df, aes(x = id,y = features.plot)) + geom_point(aes(color = avg.exp.scaled, size = pct.exp))  + scale_color_gradient(low='beige',high='tomato') + theme_classic() + theme(axis.text.x = element_text(size=18),axis.text.y=element_text(size=18),axis.title = element_blank()) + scale_size("Percent\nExpression", range = c(0.05,7)) -> p1
ggplot(AVM_vs_CTRL_art_exp_df, aes(x = id,y = features.plot)) + geom_point(aes(color = avg.exp.scaled, size = pct.exp))  + scale_color_gradient(low='beige',high='tomato') + theme_classic() + theme(axis.text.x = element_text(size=18),axis.text.y=element_blank(),axis.title = element_blank()) + scale_size("Percent\nExpression", range = c(0.05,7)) -> p2 
pdf("./Figure7h.pdf", height = 18, width = 9)
(p1|p2 ) + plot_layout(guides = 'collect')
dev.off() 

