


library(openxlsx) 
library(dplyr)
library(ggplot2) 
library(ggpmisc)

setwd("./Figure8/")



# Figure 8a
read.xlsx("./invitro_function_assay_1_2_exp.xlsx", sheet = 1) -> df1
factor(df1$group, levels = c('OE-CTRL', 'PEG10-OE-1', 'PEG10-OE-2', 'siRNA-CTRL','PEG10-siRNA-1', 'PEG10-siRNA-2')) -> df1$group 

df1$type = ifelse(grepl("OE", df1$group), "OE", 'siRNA')
df1 %>% dplyr::group_by(Day, group, type) %>% dplyr::summarise(mean_Optical_density = mean(Optical_density)) -> summary_df 
ggplot() + geom_point(data = summary_df, size = 4, alpha = 0.7,  shape = 21,aes(fill = group, x = Day %>% as.character(), y = mean_Optical_density))+ theme_bw() + theme(text = element_text(size = 20)) + ggthemes::scale_fill_tableau()+ ggthemes::scale_color_tableau() + geom_path(data = summary_df, aes(color = group, group=group,x = Day %>% as.character(), y = mean_Optical_density), size = 1.2, linetype= "dashed") + facet_wrap(~type, ncol = 1, scales = "free_y") 

df1$Compare_label <- paste0("D", df1$Day, "_", df1$group)

lapply(unique(df1$Day), function(x){
  t.test(df1$Optical_density[df1$Day %in% x & df1$group %in% c("OE-CTRL")], df1$Optical_density[df1$Day %in% x & df1$group %in% c('PEG10-OE-1','PEG10-OE-2')])$p.value -> OE_pval 
  t.test(df1$Optical_density[df1$Day %in% x & df1$group %in% c("siRNA-CTRL")], df1$Optical_density[df1$Day %in% x & df1$group %in% c('PEG10-siRNA-1','PEG10-siRNA-2')])$p.value -> si_pval 
  rbind(data.frame(type = 'OE', pval = OE_pval, Day = x), data.frame(type = "siRNA", pval = si_pval, Day = x))
}) %>% bind_rows() -> ttest_pval_stat

ttest_pval_stat$pval_label <- formatC(ttest_pval_stat$pval, digits=3)
left_join(ttest_pval_stat, summary_df %>% dplyr::group_by(type, Day) %>% dplyr::summarise(y_position = max(mean_Optical_density))) -> ttest_pval_stat 
ttest_pval_stat$adj_pval <- p.adjust(ttest_pval_stat$pval, method = 'BH')
ttest_pval_stat$adj_pval_label <- formatC(ttest_pval_stat$adj_pval, digits=3)

pdf("./Figure8a.pdf", width = 19, height = 7)
ggplot() + 
  geom_point(data = summary_df, size = 2.8, alpha = 0.9,  shape = 21,aes(fill = group, x = Day %>% as.character(), y = mean_Optical_density)) + 
  geom_path(data = summary_df, aes(color = group, group=group,x = Day %>% as.character(), y = mean_Optical_density), size = 1.2, linetype= "dashed") + 
  geom_boxplot(outlier.alpha = 0, data = df1, position = "identity", width = 0.08, size= 0.7, alpha = 0.6, aes(x = Day %>% as.character(), y = Optical_density, fill= group, group = paste0(Day, group))) + 
  geom_text(data = ttest_pval_stat, aes(x = Day %>% as.character(),y = y_position + 0.15, label = adj_pval_label), size = 5.5) +
  theme_bw() + xlab("Day") + ylab("Optical density") + 
  theme(text = element_text(size = 20), strip.text = element_text(size = 20), plot.margin = unit(c(0.6,0.6,0.6,0.6), 'cm')) + 
  facet_wrap(~type, nrow = 1, scales = "free_y") +
  ggthemes::scale_fill_tableau() + 
  ggthemes::scale_color_tableau() 
dev.off() 


# Figure 8b 
readRDS("./PEG10_Shareseq_CellCycle_meta.rds") -> shareseq_cc_meta
pdf('./Figure8b.pdf',width=5,height=6.5)
ggplot(shareseq_cc_meta, aes(x = genotype_bulk, y = S.Score + G2M.Score)) + 
  geom_violin(scale = 'width', size =0.7, aes(fill = genotype_bulk), color = "black", alpha = 0.7) + 
  geom_boxplot(outlier.alpha = 0, width = 0.25, size = 0.8) + 
  ggsignif::geom_signif(comparisons = list(c("CTRL","PEG10_RNAi"), c("CTRL", 'PEG10_OE'), c('PEG10_RNAi', "PEG10_OE")), step_increase = 0.1, test = "t.test", textsize = 7) + 
  theme_classic() + theme(text = element_text(size = 20), axis.title.x = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), 'cm'), legend.position = 'non2') + 
  ylim(c(min(shareseq_cc_meta$S.Score + shareseq_cc_meta$G2M.Score), 1.5*max(shareseq_cc_meta$S.Score + shareseq_cc_meta$G2M.Score))) + 
  scale_fill_manual(values = c("PEG10_RNAi" = "steelblue2", "CTRL" = "black", "PEG10_OE" = 'tomato')) 
dev.off() 



# Figure 8c~8f 
read.xlsx("./invitro_function_assay_1_2_exp.xlsx", sheet = 2) -> df1
factor(df1$group, levels = c('OE-CTRL', 'PEG10-OE-1', 'PEG10-OE-2', 'siRNA-CTRL','PEG10-siRNA-1', 'PEG10-siRNA-2')) -> df1$group

melt(df1, id.vars = c("group","batch")) -> df1

ggplot(df1[grepl("OE", df1$group) & df1$variable %in% 'Junction_count', ], aes(x = group, y = value)) + geom_boxplot(outlier.alpha = 0, size = 0.8) + geom_point(size = 5) + theme_bw() + theme(text = element_text(size = 20), axis.title.x = element_blank(), strip.text.x = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust= 1)) + ggsignif::geom_signif(comparisons = list(c('OE-CTRL', 'PEG10-OE-1'),c('OE-CTRL', 'PEG10-OE-2')), test = "t.test", step_increase = 0.1, textsize = 6) + facet_wrap(~variable, nrow = 1, scales = "free_y") + ylab("Number") + ylim(c(min(df1$value[grepl("OE", df1$group) & df1$variable %in% 'Junction_count']), 1.1*max(df1$value[grepl("OE", df1$group) & df1$variable %in% 'Junction_count']))) -> p1

ggplot(df1[grepl("OE", df1$group) & df1$variable %in% 'Tube_count', ], aes(x = group, y = value)) + geom_boxplot(outlier.alpha = 0, size = 0.8) + geom_point(size = 5) + theme_bw() + theme(text = element_text(size = 20), axis.title.x = element_blank(), strip.text.x = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust= 1)) + ggsignif::geom_signif(comparisons = list(c('OE-CTRL', 'PEG10-OE-1'),c('OE-CTRL', 'PEG10-OE-2')), test = "t.test", step_increase = 0.1, textsize = 6) + facet_wrap(~variable, nrow = 1, scales = "free_y")+ ylim(c(min(df1$value[grepl("OE", df1$group) & df1$variable %in% 'Tube_count']), 1.15*max(df1$value[grepl("OE", df1$group) & df1$variable %in% 'Tube_count']))) + ylab('Number') -> p2

ggplot(df1[grepl("siRNA", df1$group) & df1$variable %in% 'Junction_count', ], aes(x = group, y = value)) + geom_boxplot(outlier.alpha = 0, size = 0.8) + geom_point(size = 5) + theme_bw() + theme(text = element_text(size = 20), axis.title.x = element_blank(), strip.text.x = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust= 1)) + ggsignif::geom_signif(comparisons = list(c('siRNA-CTRL', 'PEG10-siRNA-1'),c('siRNA-CTRL', 'PEG10-siRNA-2')), test = "t.test", step_increase = 0.1, textsize = 6) + facet_wrap(~variable, nrow = 1, scales = "free_y") + ylab("Number") + ylim(c(min(df1$value[grepl("siRNA", df1$group) & df1$variable %in% 'Junction_count']), 1.1*max(df1$value[grepl("siRNA", df1$group) & df1$variable %in% 'Junction_count']))) -> p3

ggplot(df1[grepl("siRNA", df1$group) & df1$variable %in% 'Tube_count', ], aes(x = group, y = value)) + geom_boxplot(outlier.alpha = 0, size = 0.8) + geom_point(size = 5) + theme_bw() + theme(text = element_text(size = 20), axis.title.x = element_blank(), strip.text.x = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust= 1)) + ggsignif::geom_signif(comparisons = list(c('siRNA-CTRL', 'PEG10-siRNA-1'),c('siRNA-CTRL', 'PEG10-siRNA-2')), test = "t.test", step_increase = 0.1, textsize = 6) + facet_wrap(~variable, nrow = 1, scales = "free_y") + ylab("Number") + ylim(c(min(df1$value[grepl("siRNA", df1$group) & df1$variable %in% 'Tube_count']), 1.1*max(df1$value[grepl("siRNA", df1$group) & df1$variable %in% 'Tube_count']))) + ylab('Number') -> p4

pdf("./Figure8c_8f.pdf", width = 20, height = 6)
list(p1,p2,p3,p4) %>% patchwork::wrap_plots(ncol = 4, byrow = T) + theme(plot.margin = unit(c(0.8,0.8,0.8,0.8) , 'cm'))
dev.off() 


