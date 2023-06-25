
library(Seurat) 
library(ArchR) 
library(dplyr) 
library(ggplot2) 
library(ggthemes) 
library(reshape2) 
setwd("./Figure1/")

# Figure 1b 
readRDS('./scRNA_total_cell_metadata.rds') -> plc_total_cell_metadata
pdf('./Figure1b.pdf',width=7.5,height=6)
ggplot(plc_total_cell_metadata, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 0.2, aes(color = low_resolution_anno)) + theme_classic() + ggthemes::scale_color_tableau(palette = "Tableau 20")
dev.off() 

# Figure 1c 
factor(plc_total_cell_metadata$origin, levels = (c("Fetal","Maternal","Mix"))) -> plc_total_cell_metadata$origin 
plc_total_cell_metadata <- plc_total_cell_metadata %>% arrange(origin)

pdf('./Figure1c.pdf',width=10,height=9)
ggplot(plc_total_cell_metadata, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 0.2, aes(color = origin)) + theme_classic() + scale_color_manual(values = c("Fetal" = "cornflowerblue","Maternal" = "orange", "Mix" = "grey")) + theme(text = element_text(size = 20), legend.position = "right") 
dev.off()


# Figure 1d 
factor(plc_total_cell_metadata$Pheno_label, levels = c("CTRL",'PE','iTSC')) -> plc_total_cell_metadata$Pheno_label
plc_total_cell_metadata <- plc_total_cell_metadata %>% arrange(Pheno_label)

pdf('./Figure1d.pdf',width=10,height=9)
ggplot(plc_total_cell_metadata, aes(x = UMAP_1, y = UMAP_2)) + geom_point(size = 0.2, aes(color = Pheno_label)) + theme_classic() + scale_color_manual(values = c("iTSC" = "limegreen","CTRL" = "grey", "PE" = "pink")) + theme(text = element_text(size = 20), legend.position = "right")  
dev.off()


# Figure 1e 
atac_umap_color <- c('VCT/VCTp' = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[1], 
                     "SCT" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[3], 
                     "EVT" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[4], 
                     "Endo" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[5], 
                     "Fibro" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[6], 
                     "Myeloid" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[8], 
                     "CD4" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[14], 
                     "CD8" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[17], 
                     "B" = ggthemes::ggthemes_data$tableau[["color-palettes"]]$regular$`Tableau 20`$value[12])


readRDS("./ATAC_overivew_metadata.rds") -> umap_df 
pdf('./Figure1e.pdf',width=7,height=7)
ggplot(umap_df, aes(x = `IterativeLSI#UMAP_Dimension_1`, y= `IterativeLSI#UMAP_Dimension_2`, color = Integrated_annotation)) + 
  geom_point(size = 0.5, alpha = 1) + 
  scale_color_manual(values = atac_umap_color) + 
  theme_classic() + 
  xlab("") + ylab("") + 
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
dev.off()


# Figure 1f 
plotlist_5 <- readRDS('/scRNA_celltype_marker_expression_list.rds')

plotlist <- list()
for(N in seq(1:length(plotlist_5))){
  if(N < 16){
    plotlist_5[[N]] -> x
    x$labels$title -> titleid 
    x + 
      theme(axis.text.x= element_blank(), axis.title.x = element_blank(),axis.text.y=element_blank(), axis.title.y = element_text(size=20,angle=0,hjust=-5,vjust=0), 
            plot.title = element_blank(), legend.position = 'none') + 
      ylab(titleid) + ggthemes::scale_fill_tableau(palette = "Tableau 20") -> plotlist[[N]]
  } else{
    plotlist_5[[N]] -> x
    x$labels$title -> titleid 
    x + 
      theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = 0.5, size = 25), axis.title.x = element_blank(),axis.text.y=element_blank(), 
            axis.title.y = element_text(size=20,angle=0,hjust=-5,vjust=0), plot.title = element_blank(), legend.position = 'none') + 
      ylab(titleid) + ggthemes::scale_fill_tableau(palette = "Tableau 20") -> plotlist[[N]]
  }
}

pdf('./Figure1f.pdf',width=7,height=10)
patchwork::wrap_plots(plotlist = plotlist, ncol = 1) 
dev.off() 


# Figure 1h, inset Scissor 
readRDS('./Placenta_total_cell_downsample_1_scissor_meta.rds') -> scissor_downsample_meta_df1
table(scissor_downsample_meta_df1$low_resolution_anno, scissor_downsample_meta_df1$scissor)-> tb
tb/rowSums(tb) -> tb
tb %>% as.data.frame() -> scissor_plot_df 
scissor_plot_df$Celltype = scissor_plot_df$Var1
scissor_plot_df$Celltype <- as.character(scissor_plot_df$Celltype)
dcast(scissor_plot_df, Celltype~Var2, fun.aggregate = mean, value.var = "Freq") -> order_df
order_df %>% arrange(Scissor_pos, Scissor_neg) -> order_df
factor(scissor_plot_df$Celltype, levels = (order_df$Celltype)) -> scissor_plot_df$Celltype
scissor_plot_df$Var2 -> scissor_plot_df$Scissor_Type
factor(scissor_plot_df$Scissor_Type, levels = rev(c("Scissor_pos","Scissor_neg",'background_cells'))) -> scissor_plot_df$Scissor_Type
scissor_plot_df$plot_id <- ifelse(scissor_plot_df$Scissor_Type %in% "background_cells","Background", ifelse(scissor_plot_df$Scissor_Type %in% 'Scissor_pos', "PE-associated","CTRL-associated"))
factor(scissor_plot_df$plot_id, levels = c("Background","CTRL-associated","PE-associated")) -> scissor_plot_df$plot_id 

pdf("./Figure1h_Scissor.pdf", width = 6.5, height = 5)
ggplot(scissor_plot_df, aes(x = (Celltype), y = Freq)) + geom_col(color = "black", size = 0.8, alpha = 0.7, aes(fill = plot_id), position = "stack") + scale_fill_manual(values = rev(c("gold2", "cornflowerblue", "dimgrey"))) + theme_classic() + theme(text = element_text(size = 20), axis.title = element_blank(), legend.position = "right", legend.title = element_blank(), plot.margin = unit(c(0.8,0.8,0.8,0.8), "cm")) + coord_flip() + scale_y_continuous(position = "right")
dev.off()


# Figure 1h SCAVENGE
readRDS('./Figure1h_SCAVENGE_result_with_meta.rds') ->  scavenger_res
scavenger_res <- left_join(scavenger_res,meta %>% dplyr::select(Sample,Tissue_type,Phenotype))
table(scavenger_res$Integrated_annotation,scavenger_res$TRS>quantile(scavenger_res$TRS)[4]) -> tb1
(tb1/rowSums(tb1) )%>% melt -> df1
df1$Var1 <- factor(df1$Var1,levels=arrange(df1[df1$Var2,],-value)$Var1%>%as.character())
ggplot(df1,aes(x=Var1,fill=Var2,y=value)) + geom_col(size=1,color='black',alpha=0.8) + theme_classic()  + scale_fill_manual(values=c('#868686FF','#EFC000FF')) + ylab('PE GWAS Associated') + theme(text = element_text(size=20),axis.text.x = element_text(angle=90,hjust=1,vjust=0),axis.title.x = element_blank())



# Figure 1h combine SCAVENGE and Scissor 
scissor_meta_down <- readRDS("./Placenta_total_cell_downsample_2_scissor_meta.rds")

table(scissor_meta_down$low_resolution_anno,scissor_meta_down$Scissor_Type %in% 'Scissor_pos') -> tb3
(tb3/rowSums(tb3) )%>% melt -> df3

df4 <- left_join(df3,df1,by=c('Var1','Var2'))
df4$celltype_large <- 'Tropho'
df4$celltype_large[df4$Var1 %in% 'Endothelia'] <- 'Endo'
df4$celltype_large[df4$Var1 %in% 'Fibroblast'] <- 'Stromal'
df4$celltype_large[df4$Var1 %in% c('T cell','B cell','NK','M2','HB','Monocyte')] <- 'Immune'

ggplot(df4[df4$Var2 %in% TRUE,],aes(x=value.x,y=value.y)) + geom_point(aes(fill=celltype_large),pch=21,size=5) + ggrepel::geom_text_repel(mapping=aes(label=Var1),size=6) + xlab('Scissor / Phenotype / scRNA') + ylab('Scavenge / GWAS / scATAC') + theme_bw() + theme(text = element_text(size=20)) + scale_x_sqrt() + scale_y_sqrt() + scale_fill_manual(values=c('red','cornflowerblue','forestgreen','orange'))
