
library(data.table)
library(visNetwork)
library(reshape2)
library(RcisTarget)
library(dplyr)
library(ArchR)
library(pheatmap) 
library(ComplexHeatmap) 
library(openxlsx)
library(circlize)
library(patchwork)

setwd('./Fiture3/')

# Figure 3a 
data(motifAnnotations_hgnc)
motifRankings <- importRankings("./hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather") # note not using the long one
load("./selected_features_Velocity_genes.RData") # selected_features
load('./ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix.20210811.Rdata') # ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix

ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel <- ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$p_val_adj<0.05,]

t1 <- cisTarget(selected_features, motifRankings,
                motifAnnot=motifAnnotations_hgnc,verbose = T)
t1$tf_simple <- gsub(' .+','',t1$TF_highConf)

t1$NMF <- 'nothing'
selected_features[selected_features %in% t1$tf_simple] -> tf_in_selection
selected_reduced <- tf_in_selection
t1$NMF[t1$tf_simple  %in% selected_reduced] <- 'good'

signifMotifNames_list_test <- t1[t1$NMF %in% "good",] %>% dplyr::select(-enrichedGenes, -rankAtMax, -TF_lowConf, -geneSet, -nEnrGenes) %>% dplyr::filter(NES>3) %>% dplyr::select(motif) %>% unlist() 

t1[t1$NMF %in% "good"] %>% dplyr::select(-enrichedGenes, -rankAtMax, -TF_lowConf, -geneSet, -nEnrGenes) %>% dplyr::filter(NES>3) %>% dplyr::select(motif, TF_highConf) %>% as.data.frame() -> motif_TF_name_df 

gene <- selected_features

incidenceMatrix <- getSignificantGenes(gene, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames_list_test,
                                       plotCurve=FALSE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix
# 
# 
openxlsx::read.xlsx("./tropho_sub_lineage_mobile_score.xlsx") -> anno_df

anno_df$anno_genes[anno_df$anno_type %in% "VCT"] -> anno_vct
anno_df$anno_genes[anno_df$anno_type %in% "VCTp"] -> anno_vctp
anno_df$anno_genes[anno_df$anno_type %in% "SCT"] -> anno_sct
anno_df$anno_genes[anno_df$anno_type %in% "EVT"] -> anno_evt

interesting_genes <- readRDS("/Users/wanjin/Desktop/SynologyDrive/EulerianTechnology/SingleCell/Placenta/results/20210811/Tropho_velocity_interesting_gene.rds")

incidenceMatrix -> incidenceMatrix_dup
rownames(incidenceMatrix_dup) %>% as.data.frame() -> motif_name
colnames(motif_name) <- 'motif' 
left_join(motif_name,motif_TF_name_df ) -> motif_name
motif_name$TF_highConf <- gsub("\\(.+","", motif_name$TF_highConf)
motif_name$TF_highConf <- gsub("; ","/", motif_name$TF_highConf)
motif_name$TF_highConf <- gsub(" ","", motif_name$TF_highConf)
motif_name %>% as.matrix() -> motif_name_mtx
rownames(motif_name_mtx) <- motif_name_mtx[,1]
motif_name_mtx <- motif_name_mtx[,2]
cbind(motif_name_mtx, incidenceMatrix_dup) -> test
rownames(test) <- test[,1]
test <- test[,c(2:ncol(test))]
test22 <- test[rownames(test) %in% c("NFYB","EZH2","YY1","PBX1"), ]
rownames(test22) <- paste0("motif_", rownames(test22))

edges <- melt(test22)
edges <- edges[which(edges[,3]==1),1:2] ## remove genes have no linkage with any motif
colnames(edges) <- c("from","to")

motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))

nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("circle", length(genes))))

edges <- unique(edges)
nodes<-unique(nodes)

# node border width 
nodes$borderWidth <- 2.5

nodes$label <- gsub(".+_","",nodes$label)
nodes$id <- gsub(".+_","",nodes$id)
nodes$title <- gsub(".+_","",nodes$title)
nodes$shape[nodes$id %in% c("YY1","PBX1","NFYB","EZH2")] <- "diamond"
nodes <- unique(nodes)
edges$from <- gsub(".+_","",edges$from) 
unique(edges) -> edges

nodes$level <- "1"
nodes$level[nodes$id %in% c(interesting_genes, "PBX1")] <- "2"

# node border color 
nodes$color.border = 'lightgrey'
nodes$color.border[nodes$id %in% c(interesting_genes, "PBX1")] <- "#424242"
    
# node background color 
nodes$color.background <- "#E0E0E0"
nodes$color.background[nodes$id %in% anno_evt] <- "#E57373"
nodes$color.background[nodes$id %in% anno_sct] <- "#FFE082"
nodes$color.background[nodes$id %in% anno_vct] <- "#A5D6A7"
nodes$color.background[nodes$id %in% anno_vctp] <- "cornflowerblue"
nodes$color.background[grepl("NFYB", nodes$id)] <- "tomato"
nodes$color.background[grepl("PBX1", nodes$id)] <- "orange"
nodes$color.background[grepl("YY1", nodes$id)] <- "lightgrey"
nodes$color.background[grepl("EZH2", nodes$id)] <- "#90A4AE"
                    
# font color 
nodes$font.color = "#E0E0E0"
nodes$font.color[nodes$id %in% c(interesting_genes, "PBX1")] <- "black"
                      
# font size 
nodes$font.size <- 10 
nodes$font.size[nodes$id %in% c(interesting_genes)] <- 50
nodes$font.size[nodes$id %in% c("PBX1","YY1","NFYB","EZH2")] <- 95
                    
# edge width 
edges$width <- 1.5 
edges$width[edges$from == edges$to] <- 3
                    
# 
visNetwork(nodes,edges) %>% 
  visOptions(selectedBy = "level", highlightNearest = T) %>% 
  visNodes(color = list("alpha" = 0.5), mass = 3) %>% 
  visEdges(arrows='to',color = "black") %>% 
  visInteraction(dragNodes = F, dragView = F, zoomView = F) %>% 
  visLayout(randomSeed = 12)  
                    
                    
# Figure 3b 
incidenceMatrix -> incidenceMatrix_dup 
rownames(incidenceMatrix_dup) %>% as.data.frame() -> motif_name
colnames(motif_name) <- 'motif' 
left_join(motif_name,motif_TF_name_df ) -> motif_name
motif_name$TF_highConf <- gsub("\\(.+","", motif_name$TF_highConf)
motif_name$TF_highConf <- gsub("; ","/", motif_name$TF_highConf)
motif_name$TF_highConf <- gsub(" ","", motif_name$TF_highConf)
motif_name %>% as.matrix() -> motif_name_mtx
rownames(motif_name_mtx) <- motif_name_mtx[,1]
motif_name_mtx <- motif_name_mtx[,2]
cbind(motif_name_mtx, incidenceMatrix_dup) -> test
rownames(test) <- test[,1]
test <- test[,c(2:ncol(test))]
rownames(test) <- paste0("motif_", rownames(test))

edges <- melt(test)
# edges <- edges[,1:2]
edges <- edges[which(edges[,3]==1),1:2] ## remove genes have no linkage with any motif
colnames(edges) <- c("from","to")
                    
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
                    
nodes <- data.frame(id=c(motifs, genes), 
                    label=c(motifs, genes), 
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("circle", length(genes))), 
                    color=c(rep("red", length(motifs)), rep("cornflowerblue", length(genes))))
                    
edges <- unique(edges)
nodes<-unique(nodes)
                    
# Visnetwork for Motif 
nodes[nodes$shape!='circle' | nodes$id %in% tf_in_selection,] -> TF_nodes
edges[edges$from %in% TF_nodes$id & edges$to %in% TF_nodes$id,] -> TF_edges

TF_nodes$id <- gsub('.+_','',TF_nodes$id)
TF_nodes$label <- gsub('.+_','',TF_nodes$label)
TF_nodes$title <- gsub('.+_','',TF_nodes$title)
TF_edges$from <- gsub('.+_','',TF_edges$from)
TF_edges$to <- gsub('.+_','',TF_edges$to)

TF_nodes <- TF_nodes[TF_nodes$shape %in% 'circle',] %>% unique() 
TF_edges <- TF_edges %>% unique() 
# TF_nodes$color <- 'gold'
TF_edges$color <- 'gray'
TF_edges$color[TF_edges$from == TF_edges$to] <- 'red'
TF_edges$width=0.8
TF_edges$width[TF_edges$from == TF_edges$to] <- 3
TF_edges$color <- factor(TF_edges$color,levels=c('gray','red'))
TF_edges <- arrange(TF_edges,color)
    
TF_gene_paired_flow_counts <- left_join(TF_edges[TF_edges$from != TF_edges$to,] %>% dplyr::group_by(from) %>% dplyr::summarise(gene=unique(from),out_arrows=n()),TF_edges[TF_edges$from != TF_edges$to,] %>% dplyr::group_by(to) %>% dplyr::summarise(gene=unique(to),in_arrows=n()))
    
TF_gene_paired_flow_counts$from <- NULL
TF_gene_paired_flow_counts$to <- NULL
TF_gene_paired_flow_counts$flow_sum = TF_gene_paired_flow_counts$out_arrows - TF_gene_paired_flow_counts$in_arrows
    
TF_gene_paired_flow_counts$TF_hierachy_levels <- as.numeric(factor(TF_gene_paired_flow_counts$flow_sum))
    
TF_gene_paired_flow_counts$id <- TF_gene_paired_flow_counts$gene
    
TF_nodes <- left_join(TF_nodes,TF_gene_paired_flow_counts)
TF_nodes$color[TF_nodes$TF_hierachy_levels == 1] <- 'darkgray'
TF_nodes$color[TF_nodes$TF_hierachy_levels == 2] <- 'orange'
TF_nodes$color[TF_nodes$TF_hierachy_levels == 3] <- 'tomato'
TF_nodes$color.border = 'black'
TF_nodes$borderWidth <- 2.5
TF_nodes$color.background <- TF_nodes$color
TF_nodes$color <- NULL 
TF_nodes$size=40
            
TF_nodes <- arrange(TF_nodes, TF_hierachy_levels)
data.frame(x = c(-2,2,3.5,-3.5,0,0), y = c(0,0,-2,-2,-2,-4), id = c("EZH2","YY1","NFYA","NFYC","PBX1","NFYB")) -> pos_df
left_join(pos_df, TF_nodes) -> TF_nodes
pos_df$id ->  rownames(pos_df)
pos_df %>% dplyr::select(-id) %>% as.matrix() -> pos_mtx
            
visNetwork(TF_nodes,TF_edges) %>% 
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = pos_mtx) %>% 
  visNodes(font=list('color'='white','size'=20), 
           x = TF_nodes$x, 
           y = TF_nodes$y, 
           fixed = T, 
           physics = F) %>% 
  visInteraction(dragNodes = F, dragView = F, zoomView = F) %>% 
  visEdges(arrows='to',color = list('alpha'=0.5))   


# Figure 3c 

ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel <- ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$p_val_adj<0.05 & ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$avg_logFC>0.3,]
data.frame(gene=selected_features,velocity=1) -> selected_features_df 
selected_features_df <- left_join(selected_features_df,ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel %>% dplyr::select(gene,avg_logFC,cluster))

large_selected_features <- unique(c(selected_features,ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene))

# re build incidence matrix  
t1 <- cisTarget(large_selected_features, motifRankings,
                motifAnnot=motifAnnotations_hgnc,verbose = T)

t1$tf_simple <- gsub(' .+','',t1$TF_highConf)
t1$tf_simple <- gsub(';','',t1$tf_simple)

t1$NMF <- 'nothing'
selected_features[selected_features %in% t1$tf_simple] -> tf_in_selection
selected_reduced <- tf_in_selection
t1$NMF[t1$tf_simple  %in% selected_reduced] <- 'good'

signifMotifNames_list_test <- t1[t1$NMF %in% "good",] %>% dplyr::select(-enrichedGenes, -rankAtMax, -TF_lowConf, -geneSet, -nEnrGenes) %>% dplyr::filter(NES>3) %>% dplyr::select(motif) %>% unlist()

t1[t1$NMF %in% "good"] %>% dplyr::select(-enrichedGenes, -rankAtMax, -TF_lowConf, -geneSet, -nEnrGenes) %>% dplyr::filter(NES>3) %>% dplyr::select(motif, TF_highConf) %>% as.data.frame() -> motif_TF_name_df

gene <- selected_features

incidenceMatrix <- getSignificantGenes(gene,
                                       motifRankings,
                                       signifRankingNames=signifMotifNames_list_test,
                                       plotCurve=FALSE, maxRank=5000,
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

incidenceMatrix -> incidenceMatrix_dup
rownames(incidenceMatrix_dup) %>% as.data.frame() -> motif_name
colnames(motif_name) <- 'motif' 
left_join(motif_name,motif_TF_name_df ) -> motif_name
motif_name$TF_highConf <- gsub("\\(.+","", motif_name$TF_highConf)
motif_name$TF_highConf <- gsub("; ","/", motif_name$TF_highConf)
motif_name$TF_highConf <- gsub(" ","", motif_name$TF_highConf)
motif_name %>% as.matrix() -> motif_name_mtx
rownames(motif_name_mtx) <- motif_name_mtx[,1]
motif_name_mtx <- motif_name_mtx[,2]
cbind(motif_name_mtx, incidenceMatrix_dup) -> test
rownames(test) <- test[,1]
test <- test[,c(2:ncol(test))]
test22 <- test[rownames(test) %in% c("NFYB","EZH2","YY1","PBX1",'NFYA','NFYC'), ]
rownames(test22) <- paste0("motif_", rownames(test22))

# stats 
edges_stat <- melt(test22)
edges_stat <- edges_stat[which(edges_stat[,3]==1),1:2] ## remove genes have no linkage with any motif
colnames(edges_stat) <- c("from","to")
edges_stat$width <- 1

edges_stat <- unique(edges_stat)

edges_stat$cell_type <- ''
edges_stat$cell_type[edges_stat$to %in% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'VCTp'],anno_vctp),'PCNA','AURKA','AURKB','CENPC','PLK2','PLK1','HMGB2')] <- 'VCTp'
edges_stat$cell_type[edges_stat$to %in% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'VCT'],anno_vct)) & edges_stat$to %ni% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'VCTp'],anno_vctp),'PCNA','AURKA','AURKB','CENPC','PLK2','PLK1','HMGB2')] <- 'VCT'
edges_stat$cell_type[edges_stat$to %in% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'EVT'],anno_evt)) & edges_stat$to %ni% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'VCTp'],anno_vctp),'PCNA','AURKA','AURKB','CENPC','PLK2','PLK1','HMGB2')] <- 'EVT'
edges_stat$cell_type[edges_stat$to %in% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'SCT'],anno_sct)) & edges_stat$to %ni% na.omit(c(ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix_sel$gene[ctrl_trophoblast_four_major_classes_velocity_gene_expression_matrix$cluster %in% 'VCTp'],anno_vctp),'PCNA','AURKA','AURKB','CENPC','PLK2','PLK1','HMGB2')] <- 'SCT'
edges_stat$from <- gsub('motif_','',edges_stat$from)
edges_stat <- unique(edges_stat)

edges_stat %>% dplyr::group_by(to,cell_type) %>% dplyr::summarise('TF'=paste0(sort(unique(from)),collapse=";")) -> tt1
tt1$max_TF_layer <- ''
tt1$max_TF_layer[!grepl('EZH2|YY1',tt1$TF) & grepl('NFYA|NFYC|PBX1',tt1$TF)] <- '2 Medium: NFYA/C;PBX1'
tt1$max_TF_layer[grepl('EZH2|YY1',tt1$TF) & !grepl("NFY|PBX",tt1$TF)] <- '3 Late: PRC2: EZH2/YY1'
tt1$max_TF_layer[grepl('NFYB',tt1$TF) ] <- '1 Early: NFYB'
table(tt1$cell_type,tt1$max_TF_layer) -> tb3
tb3 <- tb3[2:nrow(tb3),]
tb3 <- melt(tb3) %>% dcast(Var1~Var2,fun.aggregate=sum,value.var='value')
rownames(tb3) <- tb3$Var1
tb3 <- tb3[,c(2:ncol(tb3))]
tb3 <- tb3[,colSums(tb3)>0]
# t(tb3) -> tb3
# tb3 <- t(tb3)
tb3 <- tb3/rowSums(tb3)
tb3 <- t(tb3)

library(ComplexHeatmap)
tb3[,c('VCTp','VCT','SCT','EVT')] -> tb4
# scale by row 
tb4 <- t(tb4)
scale(tb4) -> tb4 
tb4 -> tb5
tb5[tb5>1] <- 1
tb5[tb5< -1] <- -1

pheatmap::pheatmap(t(tb5),cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c('blueviolet','black','gold'))(30),border_color = 'white',cellwidth = 40,cellheight = 40)


# Figure 3d 
z_score_mtx_list <- readRDS("./ATAC_Z_score_hm_mtx_list.rds")

col_order <- c('VCTp_4','VCTp_8','VCTp_1','VCT_2', 'VCT_5','VCT_3','VCT_7','VCT_10','VCT_0','SCT_6','SCT_12','EVT_11','EVT_13','EVT_9')

genes_for_annotation <- unique(c('BACH1','BACH2','FOS_31','FOS..JUND_341','FOS..JUN_326','STAT1_49',"HEY1_244","PBX1_8","NFYB_601","SOX2_617",'GCM1_80', 'CDX1_487','HES1_497','POU5F1_310','CTCF_21','YY1_52','HOXA13_499','NFYA_305','NFYC_539','EGR3_162','EGR1_563','SNAI1_443','ZEB1_320','NR2F2_307','SPIB_619','TFAP2C_235','MEF2A_593','TEAD4_625','BCL6_486','RELA_17','NFKB1_198','BCL6_486','FOXA2_573','TCF7_623','IRF8_89'))

statistics_mtx_for_TF_in_each_celltype_for_plot_2 <- z_score_mtx_list[['Z_Score']]
statistics_mtx_for_TF_in_each_celltype_for_plot_1 <- z_score_mtx_list[['Log.Adj.Pval']]


hclust(dist(statistics_mtx_for_TF_in_each_celltype_for_plot_2),method = 'ward.D2') -> clust_obj
hclust(dist(t(statistics_mtx_for_TF_in_each_celltype_for_plot_2)),method = 'ward.D2') -> clust_obj2
cutree(clust_obj,5)[clust_obj$labels[clust_obj$order]] -> split_row_obj

col_order %>% as.data.frame -> split_col_df 
colnames(split_col_df) <- "cellgroup"
split_col_df$group <- "2" 
split_col_df$group[grepl("VCTp", split_col_df$cellgroup)] <- "1" 
split_col_df$group[grepl("SCT", split_col_df$cellgroup)] <- "3" 
split_col_df$group[grepl("EVT", split_col_df$cellgroup)] <- "3" 
split_col_df$cellgroup -> rownames(split_col_df) 
split_col_df %>% dplyr::select(group) -> split_col_df

statistics_mtx_for_TF_in_each_celltype_for_plot_2[clust_obj$labels[clust_obj$order],clust_obj2$labels[clust_obj2$order]] -> statistics_mtx_for_TF_in_each_celltype_for_plot_2
statistics_mtx_for_TF_in_each_celltype_for_plot_1[clust_obj$labels[clust_obj$order],clust_obj2$labels[clust_obj2$order]] -> statistics_mtx_for_TF_in_each_celltype_for_plot_1

right_annotation_link = rowAnnotation(link = anno_mark(at = which(grepl(paste0("z:",genes_for_annotation,collapse = "|"),rownames(statistics_mtx_for_TF_in_each_celltype_for_plot_2))), 
                                                       labels = gsub('z:','',gsub('_.+','',rownames(statistics_mtx_for_TF_in_each_celltype_for_plot_2)[grepl(paste0("z:",genes_for_annotation,collapse = "|"),rownames(statistics_mtx_for_TF_in_each_celltype_for_plot_2))])), 
                                                       labels_gp = gpar(fontsize = 20), padding = unit(1, "mm")))

col_fun_Z = colorRamp2(c(-1,-0.8,-0.3, 0, 0.3, 0.8,1), c("darkorchid4","darkorchid1","black","black","black","gold4",'gold'))
col_fun_log = colorRamp2(c(0,0.8,1,1.3,3), c("black","black","black",'gold1','gold1'))
top_annotation_block <- HeatmapAnnotation(CellType = anno_block(gp = gpar(fill = c('white')), labels = c('VCTp','VCT','SCT/EVT'),labels_gp = gpar(col = "black", fontsize = 25)))

Heatmap(
  statistics_mtx_for_TF_in_each_celltype_for_plot_2[clust_obj$labels[clust_obj$order],col_order],
  cluster_rows=F, cluster_columns=F, show_row_names=FALSE, show_column_names=T, 
  col = col_fun_Z, name = 'Z', height = unit(20,"cm"), width = unit(15,"cm"), column_names_gp = gpar(fontsize = 25), 
  column_title = "PE over CTRL Z", column_title_side = 'top',  column_title_gp = gpar(fontsize = 30), 
  split = split_row_obj,
  column_split = split_col_df,
  right_annotation = right_annotation_link, 
  top_annotation = top_annotation_block
) -> ph_z 

Heatmap(
  statistics_mtx_for_TF_in_each_celltype_for_plot_1[clust_obj$labels[clust_obj$order],col_order],
  cluster_rows=F, cluster_columns=F, show_row_names=FALSE, show_column_names=T, 
  col = col_fun_log, name = 'log.adj.Pval', height = unit(20,"cm"), width = unit(15,"cm"), column_names_gp = gpar(fontsize = 25), 
  column_title = "PE over CTRL Log adj Pval",  column_title_side = 'top', column_title_gp = gpar(fontsize = 30), 
  split = split_row_obj,
  column_split = split_col_df, 
  top_annotation = top_annotation_block
) -> ph_log

pdf(paste0('./Figure3d.pdf'),height=15,width=15)
print(ph_log+ph_z)
dev.off() 

# Figure 3e 
tf_bulk_atac_hint_analysis_result_df <- readRDS("Combine_TF_bulk_atac_hint_result.rds") 
tf_differential_analysis_res <- read.delim('./differential_statistics.txt')  
tf_differential_analysis_res$TC_CTRL <- as.numeric(as.character(tf_differential_analysis_res$TC_CTRL)) 
tf_differential_analysis_res$TC_PE <- as.numeric(as.character(tf_differential_analysis_res$TC_PE))
tf_differential_analysis_res$P_values <- as.numeric(as.character(tf_differential_analysis_res$P_values)) 
tf_differential_analysis_res$Z_score <- as.numeric(as.character(tf_differential_analysis_res$Z_score))

tf_differential_analysis_res$protection_ratio <- tf_differential_analysis_res$TC_CTRL/tf_differential_analysis_res$TC_PE
tf_differential_analysis_res$TF_name <- gsub('.+\\.','',gsub('\\(.+','',tf_differential_analysis_res$Motif))

ggplot(tf_differential_analysis_res,aes(x=log(protection_ratio),y=-log(P_values))) + geom_point(aes(),color='gray',size=4) + geom_point(aes(color=abs(Z_score)>3),size=3) + theme_classic() + scale_color_manual(values=c('black','red')) + ggrepel::geom_text_repel(data=tf_differential_analysis_res[abs(tf_differential_analysis_res$Z_score)>3,],segment.alpha = 0,aes(label=TF_name),min.segment.length = unit(3,'cm'),size=8) + xlim(c(-0.2,0.2)) + theme(legend.position = 'bottom',text = element_text(size=20)) -> p1

tf_motif_name <- unique(tf_bulk_atac_hint_analysis_result_df$file[grepl('NFY|Dux|ZNF384',tf_bulk_atac_hint_analysis_result_df$file)])
data.frame(tf_motif_name=tf_motif_name,tf_name=gsub('.+\\.','',gsub('.txt','',tf_motif_name))) -> tf_motif_name_df
tf_motif_name_df$tf_name <- factor(tf_motif_name_df$tf_name,levels=c('NFYA','NFYB','Dux','ZNF384'))
tf_motif_name_df <- arrange(tf_motif_name_df,tf_name)
tf_motif_name_df$tf_motif_name <- factor(tf_motif_name_df$tf_motif_name,levels = tf_motif_name_df$tf_motif_name,ordered = T)

lapply(tf_motif_name_df$tf_name,function(x){
  tf_motif_name_id <- tf_motif_name_df$tf_motif_name[tf_motif_name_df$tf_name%in%x]
  tf_footprint_df <- tf_bulk_atac_hint_analysis_result_df[tf_bulk_atac_hint_analysis_result_df$file %in% tf_motif_name_id,]
  tf_footprint_df <- melt(tf_footprint_df,id.vars=c('pos','file'))
  ggplot(tf_footprint_df,aes(x=pos,y=value,color=variable)) + geom_line() + scale_color_manual(values=c('PE'='red','CTRL'='black')) + theme_classic() + xlim(c(-99,99)) + theme(axis.title.x = element_blank(),text=element_text(size=20)) + ylab('Footprint') + ggtitle(x)
}) %>% wrap_plots(ncol=2,guides = 'collect') + plot_annotation(theme=theme(legend.position = 'bottom',text = element_text(size=20)))  -> p2

pdf("Figure3e.pdf", width = 12, height = 6)
p1|p2  
dev.off() 

