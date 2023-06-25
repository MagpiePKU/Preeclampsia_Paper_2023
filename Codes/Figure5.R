
library(ggplot2) 
library(ArchR)
library(ggsignif) 
library(patchwork)
setwd("./Figure5/") 

# Figure 5a 
load('./Tropho_total_FindAllMarkers_after_integration_20220525.RData') # tropho_total_marker_after_integration
tropho_total_marker_after_integration[tropho_total_marker_after_integration$cluster %in% "PE" & !grepl("\\.|^MT|^RP|^HSP", tropho_total_marker_after_integration$gene), ] -> tropho_total_marker_after_integration_clean
tropho_label <- c("PEG10","PAPPA2","PAPPA","FLT1","IGFBP1","SIGLEC6","DUSP1","CCND1","GADD45A","SLC13A3","CGA","CCL4","KISS1","HOPX","HSD17B1","PGF","MALAT1","CSF3R","FOLR1","GDF15","MMP2","HLA-G","ANXA3","H2AZ1")

pdf('./Figure5a.pdf',width=8,height=8)
ggplot() + 
  geom_vline(xintercept = 0, size = 1.5, color = "dimgrey", linetype = "dashed") + 
  geom_point(data = tropho_total_marker_after_integration_clean[-log10(tropho_total_marker_after_integration_clean$p_val) < 250, ], 
             mapping = aes(x = (avg_logFC/abs(avg_logFC))*sqrt(abs(avg_logFC)), y = -log10(p_val_adj)), 
             size = 2, fill = "lightgrey", alpha = 0.5, shape = 21, colour = "black") + 
  geom_point(data = tropho_total_marker_after_integration_clean[(-log10(tropho_total_marker_after_integration_clean$p_val)) < 250 & tropho_total_marker_after_integration_clean$p_val < 0.01 & abs(tropho_total_marker_after_integration_clean$avg_logFC) > 0.2, ], 
             mapping = aes(x = (avg_logFC/abs(avg_logFC))*sqrt(abs(avg_logFC)), y = -log10(p_val_adj), fill = avg_logFC > 0), 
             size = 2, alpha = 0.7,shape = 21, colour = "black") + 
  ggrepel::geom_label_repel(data = tropho_total_marker_after_integration_clean[tropho_total_marker_after_integration_clean$p_val < 0.01 & (tropho_total_marker_after_integration_clean$avg_logFC) < -0.2 & tropho_total_marker_after_integration_clean$gene %in% tropho_label, ], 
                            mapping = aes(x = (avg_logFC/abs(avg_logFC))*sqrt(abs(avg_logFC)), y = -log10(p_val_adj), label = gene), 
                            max.time = 100000, nudge_x = -0.02, nudge_y = 0.02, min.segment.length = unit(0.05, "mm"), size = 5) + 
  ggrepel::geom_label_repel(data = tropho_total_marker_after_integration_clean[tropho_total_marker_after_integration_clean$p_val < 0.01 & (tropho_total_marker_after_integration_clean$avg_logFC) > 0.2 & tropho_total_marker_after_integration_clean$gene %in% tropho_label, ], 
                            mapping = aes(x = (avg_logFC/abs(avg_logFC))*sqrt(abs(avg_logFC)), y = -log10(p_val_adj), label = gene), 
                            max.time = 100000, nudge_x = 0.02, nudge_y = 0.02, min.segment.length = unit(0.05, "mm"), size = 5) +
  theme_classic() + scale_fill_manual(values = c('cornflowerblue', "orange")) +
  scale_x_continuous(limits = c(-2,2), breaks = c(-1,0, 1),labels = c("PE Attenuated","0","PE Enhanced")) +
  theme(legend.position = "none", axis.title.x = element_blank(), text = element_text(size = 20), axis.ticks.x = element_blank())
dev.off()

# Figure 5b
plot_df <- readRDS("./Tropho_PEG10_expression.rds")
ggplot(plot_df[plot_df$PE_detail %ni% "DCDA_PE", ], aes(x = modified_celltype, y = PEG10)) + 
  geom_violin(scale = "width", color = "black", size = 0.9, aes(fill = modified_celltype), alpha = 0.8) + 
  geom_boxplot(outlier.alpha = 0, width = 0.15, color = "black", fill = "lightgrey", alpha = 1, size = 1) + 
  facet_wrap(~new_PE, scales = "free") + theme_classic() + 
  theme(axis.text.x= element_blank() , axis.title.x = element_blank(), legend.position = 'right', text = element_text(size = 25), legend.title = element_blank()) + scale_fill_manual(values = c("cornflowerblue","forestgreen","orange","red")) + ggtitle("CTRL vs PE Placenta") -> pp1

ggplot(plot_df[plot_df$PE_detail %in% "DCDA_PE", ], aes(x = modified_celltype, y = PEG10)) + 
  geom_violin(scale = "width", color = "black", size = 0.9, aes(fill = modified_celltype), alpha = 0.8) + 
  geom_boxplot(outlier.alpha = 0, width = 0.15, color = "black", fill = "lightgrey", alpha = 1, size = 1) + 
  facet_wrap(~tissue_type, scales = "free") + theme_classic() + 
  theme(axis.text.x=  element_blank() , axis.title.x = element_blank(), legend.position = 'right', text = element_text(size = 25), legend.title = element_blank()) + scale_fill_manual(values = c("cornflowerblue","forestgreen","orange","red")) + ggtitle("DCDA PE Placenta") -> pp2

pdf('./Figure5b.pdf',width=9,height=6)
pp1/pp2
dev.off()


# Figure 5c
ggplot(plot_df[plot_df$PE_detail %ni% "DCDA_PE" & plot_df$new_cluster %in% 'SCT_12', ], aes(x = new_PE, y = PEG10)) + geom_violin(scale = "width", color = "black", size = 0.9, aes(fill = new_PE), alpha = 0.7) + 
  geom_boxplot(outlier.alpha = 0, width = 0.15, color = "black", fill = "grey", alpha = 1, size = 1) + 
  facet_wrap(~new_cluster, scales = "free", nrow = 1) + theme_classic() + 
  theme(axis.text.x= element_text(angle = 45, hjust = 1, vjust = 1) , axis.title.x = element_blank(), legend.position = 'none', text = element_text(size = 20), legend.title = element_blank(), legend.box.just = "left") + scale_fill_manual(values = c("dimgrey","red")) + ggtitle("CTRL vs PE") + 
  geom_signif(comparisons = list(c("PE","CTRL")), test = "t.test", map_signif_level = T, tip_length = 0, size = 1, vjust = 1.5, textsize  = 8) -> pp3

ggplot(plot_df[plot_df$PE_detail %in% "DCDA_PE" &  plot_df$new_cluster %in% 'SCT_12', ], aes(x = tissue_type, y = PEG10)) + 
  geom_violin(scale = "width", color = "black", size = 0.9, aes(fill = tissue_type), alpha = 0.7) + 
  geom_boxplot(outlier.alpha = 0, width = 0.15, color = "black", fill = "lightgrey", alpha = 1, size = 1) + 
  facet_wrap(~new_cluster, scales = "free", nrow = 1) + theme_classic() + 
  theme(axis.text.x=  element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank(), legend.position = 'none', text = element_text(size = 20), legend.title = element_blank(), legend.box.just = "left") + scale_fill_manual(values = c("lightcoral","chocolate")) + ggtitle("DCDA PE") + 
  geom_signif(comparisons = list(c("normal_side","FGR_side")), test = "t.test", map_signif_level = T, tip_length = 0, size = 1, vjust = 1.5, textsize = 8) -> pp4

pdf('./Figure5c.pdf',width=4,height=8)
pp3/pp4
dev.off()

# Figure 5g 

load('./PEG10_combine_plot_df_20220421.RData')
my_comarisons=list(c('Endo_CTRL','Endo_SPE'),c('Tropho_CTRL','Tropho_SPE'),c('Endo_CTRL','Tropho_CTRL'),c('Endo_SPE','Tropho_SPE')) 
pdf("./Figure5g.pdf", width = 8, height = 8)
ggplot(PEG10_combine_plot_df,aes(x=Cell_Phenotype,y=normb)) + geom_violin(aes(fill=as.character(Cell_Phenotype)),scale='width') + scale_y_log10() + geom_boxplot(aes(fill=Phenotype),width=0.3,outlier.alpha = 0) +  scale_fill_manual(values=c('gray','darkolivegreen3','forestgreen','black','cornflowerblue','cadetblue4')) + theme_classic() + ggpubr::stat_compare_means(comparisons = my_comarisons,method = 'wilcox.test') + ylab("log PEG10 Intensity") + ggtitle('PEG10 IF')  
dev.off() 


# Figure 5h 
PEG10_dcda_df <- readRDS('./PEG10_IF_plot_for_dcda_in_figure.rds')
my_comarisons=list(c('Endo_L','Endo_S'),c('Tropho_L','Tropho_S'),c('Endo_L','Tropho_L'),c('Endo_S','Tropho_S'))

pdf("./Figure5h.pdf", width = 8, height = 8)
ggplot(PEG10_dcda_df,aes(x=Cell_Phenotype,y=normb)) + geom_violin(aes(fill=as.character(Cell_Phenotype)),scale='width') + scale_y_log10() + geom_boxplot(aes(fill=Phenotype_Pseudo),width=0.3,outlier.alpha = 0) +  scale_fill_manual(values=c('darkolivegreen3','forestgreen','gray','black','cornflowerblue','cadetblue4')) + theme_classic() + ggpubr::stat_compare_means(comparisons = my_comarisons,method = 'wilcox.test') + ylab("log PEG10 Intensity") + ggtitle('DCDA PEG10 IF')

dev.off()


# Figure 5i 
cutoff_endo = 0.1
cutoff_tropho = 0.37

PEG10_combine_plot_df_for_cutoff <- PEG10_combine_plot_df
PEG10_combine_plot_df_for_cutoff$cutoff <- 0
PEG10_combine_plot_df_for_cutoff$cutoff[PEG10_combine_plot_df_for_cutoff$Mask_type %in% 'Endo'] <- cutoff_endo
PEG10_combine_plot_df_for_cutoff$cutoff[PEG10_combine_plot_df_for_cutoff$Mask_type %in% 'Tropho'] <- cutoff_tropho
PEG10_combine_plot_df_for_cutoff$PEG10type <- ifelse(PEG10_combine_plot_df_for_cutoff$normb > PEG10_combine_plot_df_for_cutoff$cutoff,'PEG10+','PEG10-')
PEG10_combine_plot_df_for_cutoff %>% dplyr::group_by(sid,Mask_type,PEG10type) %>% dplyr::summarise(n=n()) -> PEG10_combine_plot_df_for_cutoff_summary
PEG10_combine_plot_df_for_cutoff_summary <- left_join(PEG10_combine_plot_df_for_cutoff_summary,PEG10_combine_plot_df_for_cutoff_summary %>% dplyr::group_by(sid,Mask_type) %>% dplyr::summarise(sumN=sum(n)))
PEG10_combine_plot_df_for_cutoff_summary$freq <- PEG10_combine_plot_df_for_cutoff_summary$n/PEG10_combine_plot_df_for_cutoff_summary$sumN
PEG10_combine_plot_df_for_cutoff_summary$phenotype <- ifelse(grepl('SF|PE',PEG10_combine_plot_df_for_cutoff_summary$sid),'PE','CTRL')
PEG10_combine_plot_df_for_cutoff_summary$Cell_Phenotype <- paste0(PEG10_combine_plot_df_for_cutoff_summary$phenotype,'_',PEG10_combine_plot_df_for_cutoff_summary$Mask_type) %>% factor(levels=c('CTRL_Tropho','PE_Tropho','CTRL_Endo','PE_Endo'))

pdf("./Figure5i.pdf", width = 8, height = 8)
ggplot(PEG10_combine_plot_df_for_cutoff_summary[PEG10_combine_plot_df_for_cutoff_summary$PEG10type %in% 'PEG10+' ,],aes(x=Cell_Phenotype,y=freq,fill=Cell_Phenotype))  + geom_boxplot()  +  ggbeeswarm::geom_beeswarm(aes(shape=phenotype),size=4) + theme_classic() +  scale_fill_manual(values=c('cornflowerblue','cadetblue4','darkolivegreen3','forestgreen')) + theme(text=element_text(size=12)) + ggtitle('PEG10+ Fraction') + ylab('Frequency') 
dev.off() 



