
library(dplyr) 
library(ggplot2) 
library(ArchR) 
library(ggrepel) 
library(patchwork)

setwd("./Figure6/")

# Figure 6d
arterial_PEvsCTRL_marker_df <- readRDS('./Arterial_PEvsCTRL_marker_20220523.rds')

arterial_PEvsCTRL_marker_df[!grepl("^MT|^RP|\\.", arterial_PEvsCTRL_marker_df$gene),] -> arterial_PEvsCTRL_marker_df_clean 
arterial_PEvsCTRL_marker_df_clean[arterial_PEvsCTRL_marker_df_clean$cluster %in% "PE", ] -> arterial_volcano_plot

real_cargo_gene = c("AOC1","KRT8","PAPPA2","NOTUM","KRT19","HPGD","CSH1","KRT7","TAC3","KRT18","HTRA1", "FBLN1","EBI3","MFAP5","ISM2","PRG2","PTPRF","S100P","LVRN", "PHLDA2","CST6","COL17A1", "LGALS3BP", "BGN", "ASCL2","SERPINE1", "TIMP2","PEG10")
c("PAPPA2","PEG10", "NOTCH1","NOTCH4","DLL4","HES4","JAG1","HEY1","CTNNB1","JAG2","HES1","CD36","HTRA1","FLT1","NOTUM","PAPPA","IGFBP1","MYC","CTNNB1","EGFR","PGF","TNFSF10","TIMP3","TIMP2","SAT1","IGFBP4","FN1","ITGB1","CCL4","CCL4L2","CXCL2","CCL3") -> label_genes11
unique(c(label_genes11, real_cargo_gene)) -> label_genes

pdf("./Figure6d.pdf", width = 10, height = 10)
ggplot() + 
  geom_vline(xintercept = 0, size = 1.5, linetype = "dashed", color = "grey") + 
  geom_point(data = arterial_volcano_plot, mapping = aes(x = avg_logFC, y = -log10(p_val_adj)),color = "black", size = 2.5, alpha = 0.3, shape = 21, color = "black", fill = "dimgrey") + 
  geom_point(data = arterial_volcano_plot[arterial_volcano_plot$p_val < 0.01 & (arterial_volcano_plot$avg_logFC) > 0.1, ], 
             mapping = aes(x = avg_logFC, y = -log10(p_val_adj)), shape = 21, color = "black", size = 2.5, alpha = 0.7, fill = "orange") + 
  geom_point(data = arterial_volcano_plot[arterial_volcano_plot$p_val < 0.01 & (arterial_volcano_plot$avg_logFC) < -0.1, ], 
             mapping = aes(x = avg_logFC, y = -log10(p_val_adj)), shape = 21, color = "black", size = 2.5, alpha = 0.7, fill = "cornflowerblue") + 
  ggrepel::geom_label_repel(data = arterial_volcano_plot[arterial_volcano_plot$gene %in% label_genes & arterial_volcano_plot$avg_logFC > 0 , ], 
                            mapping = aes(x = avg_logFC, y = -log10(p_val_adj), label = gene, fill = gene %in% real_cargo_gene), 
                            min.segment.length = unit(0.6, "mm"), nudge_x = 0.2, nudge_y = 0.8,size = 6) + 
  ggrepel::geom_label_repel(data = arterial_volcano_plot[arterial_volcano_plot$gene %in% label_genes & arterial_volcano_plot$avg_logFC < 0, ], 
                            mapping = aes(x = avg_logFC, y = -log10(p_val_adj), label = gene, fill = gene %in% real_cargo_gene), 
                            min.segment.length = unit(0.6, "mm"), nudge_x = -0.7, nudge_y = 0.4,size = 6) + 
  theme_classic() + scale_fill_manual(values = c("white", "#FCE4EC")) + 
  theme(legend.position = "none", text = element_text(size = 25), axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + 
  scale_x_continuous(limits = c(-2.5,2.5),breaks = c(-1.5,0,1.5), labels = c("PE Attenuated", "0", "PE Enhanced"))
dev.off() 

# Figure 6e 
plotlist_5 <- readRDS("./cargo_gene_expression_list_in_endo.rds")

color_code = c('fetal_CTRL'='cornflowerblue','fetal_PE'='navy','maternal_CTRL'='gold1','maternal_PE'='orange')

lapply(seq(1:length(plotlist_5)), function(x){
  if(x == length(plotlist_5)){
    plotlist_5[[x]]$data -> dat
    colnames(dat)[1] -> title 
    colnames(dat)[1] <- "gene" 
    dat$origin <- dat$ident
    dat$origin <- gsub("_CTRL|_PE","",dat$origin) 
    dat$PE <- dat$ident
    dat$PE <- gsub("fetal_|maternal_","",dat$PE) 
    ggplot(dat, aes(x = ident, y = gene, fill = ident)) + 
      geom_violin(scale = "width", alpha = 0.7) + 
      geom_boxplot(width = 0.1, color = "black", fill = "grey", outlier.alpha = 0) +
      theme_classic() + 
      theme(axis.text.x= element_text(angle = 45, hjust = 1, vjust = 1, size = 15), 
            axis.title.x = element_blank(), 
            axis.text.y=element_blank(), 
            axis.title.y = element_text(size=15,angle=0,hjust=-5,vjust=0), 
            plot.title = element_blank(), 
            legend.position = 'none', 
            text = element_text(size=15)) + 
      ylab(title) + 
      scale_fill_manual(values=color_code)
  }else{
    plotlist_5[[x]]$data -> dat
    colnames(dat)[1] -> title 
    colnames(dat)[1] <- "gene" 
    dat$origin <- dat$ident
    dat$origin <- gsub("_CTRL|_PE","",dat$origin) 
    dat$PE <- dat$ident
    dat$PE <- gsub("fetal_|maternal_","",dat$PE) 
    ggplot(dat, aes(x = ident, y = gene, fill = ident)) + 
      geom_violin(scale = "width", alpha = 0.7) + 
      geom_boxplot(width = 0.1, color = "black", fill = "grey", outlier.alpha = 0) +
      theme_classic() + 
      theme(axis.text.x= element_blank(), 
            axis.title.x = element_blank(), 
            axis.text.y=element_blank(), 
            axis.title.y = element_text(size=15,angle=0,hjust=-5,vjust=0), 
            plot.title = element_blank(), 
            legend.position = 'none', 
            text = element_text(size=15)) + 
      ylab(title) + 
      scale_fill_manual(values=color_code)
  }
  
}) %>% patchwork::wrap_plots(ncol =1) -> Plot_B

pdf('./Figure6e.pdf',width=4,height=6)
print(Plot_B)
dev.off()



