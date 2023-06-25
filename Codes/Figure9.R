

library(ggplot2)
library(dplyr)
library(openxlsx) 
setwd("./Figure9/")

# Figure 9a 
readRDS("./PEG10_Shareseq_DE_gene_expression.rds") -> plot_data

factor(plot_data$id, levels = rev(c("PEG10_RNAi","CTRL","PEG10_OE"))) -> plot_data$id 
pdf('./Figure9a.pdf',width=10,height= 4.5)
ggplot(plot_data,aes(x=features.plot,y = id,size = pct.exp, fill = avg.exp.scaled.mod))+
  geom_point(shape = 21, color = "dimgrey") + 
  scale_size("Percent\nExpression", range = c(2.5,10)) + 
  scale_fill_gradientn(colours = c("steelblue","white","red"),
                       guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                       name = "Average\nExpression") + 
  theme_bw() + 
  theme(text = element_text(size = 25), 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right", 
        # legend.spacing.x = unit(0.2, 'cm'),
        # legend.spacing.y = unit(0.5, 'cm'),
        legend.text.align = 0, 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20), 
        legend.key.size = unit(1, "cm"), 
        # legend.key.height = unit(0.5, "cm"), 
        legend.key.width = unit(0.7, "cm"),
        legend.direction = "vertical", plot.margin = unit(c(0.6,0.6,0.6,0.6), "cm")) 
dev.off() 
# Figure 9b 
gsea_res <- readRDS("./Significant_GSEA_result_C2_H.rds")

pdf('./Figure9b.pdf',height=7,width=13)
ggplot(gsea_res,aes(x=NES,y=pathway,fill= -sign(NES)*log10(pval))) + 
  geom_vline(xintercept = 0, size = 0.6, alpha = 0.7, color = "black") + 
  geom_col(size=0.5, alpha = 0.6,color='black')  + 
  theme_classic() + xlab('NES') + 
  theme(text=element_text(size=20),axis.title.y=element_blank(), legend.position = "right", legend.direction = 'vertical', plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), legend.text.align = 1, legend.key.size = unit(1, "cm"), legend.key.width = unit(0.7, "cm")) + 
  ggtitle("PEG10 OE GSEA", "PEG10 OE vs PEG10 RNAi") + scale_fill_gradientn(colours = c("steelblue","white","red"), guide = guide_colorbar(ticks.colour = "black",frame.colour = "black")) 
dev.off() 

# Figure 9c~9d 

read.xlsx("./invitro_function_assay_3_exp.xlsx", sheet = 1) -> df1

melt(df1, id.vars = c("group","batch")) -> df1
factor(df1$group, levels = (c('OE-CTRL','OE-CTRL+TGFb1', 'PEG10-OE-1','PEG10-OE-1+TGFb1', 'PEG10-OE-2','PEG10-OE-2+TGFb1'))) -> df1$group

df1_plot <- df1[df1$group %ni% "OE-CTRL+TGFb1", ] 

ggplot(df1_plot[df1_plot$variable %in% 'Junction_count', ], aes(x = group, y = value)) + 
  geom_boxplot(outlier.alpha = 0, size = 0.8) + 
  geom_point(size = 5) + 
  theme_bw() + xlab("Number") +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1), axis.title.x = element_blank(), strip.text = element_text(size = 20), plot.margin = unit(c(0.6,0.6,0.6,0.6), 'cm')) + 
  ggsignif::geom_signif(comparisons = list(c('PEG10-OE-1', 'PEG10-OE-1+TGFb1'),c('PEG10-OE-2', 'PEG10-OE-2+TGFb1')), test = "t.test", y_position = 260, textsize = 6, orientation = "x") + 
  ggsignif::geom_signif(comparisons = list(c('OE-CTRL', 'PEG10-OE-1'), c('OE-CTRL', 'PEG10-OE-1+TGFb1'), c('OE-CTRL', 'PEG10-OE-2'), c('OE-CTRL', 'PEG10-OE-2+TGFb1')), test = "t.test", step_increase = 0.1, textsize = 6, orientation = "x", y_position = 300) + 
  facet_wrap(~variable, nrow = 1, scales = "free_y") + ylab("Number") -> p5


ggplot(df1_plot[df1_plot$variable %in% 'Tube_count', ], aes(x = group, y = value)) + 
  geom_boxplot(outlier.alpha = 0, size = 0.8) + 
  geom_point(size = 5) + 
  theme_bw()+ xlab("Number") + 
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1), axis.title.x= element_blank(), strip.text = element_text(size = 20), plot.margin = unit(c(0.6,0.6,0.6,0.6), 'cm')) +  
  ggsignif::geom_signif(comparisons = list(c('PEG10-OE-1', 'PEG10-OE-1+TGFb1'),c('PEG10-OE-2', 'PEG10-OE-2+TGFb1')), test = "t.test", y_position = 66, textsize = 6, orientation = "x") + 
  ggsignif::geom_signif(comparisons = list(c('OE-CTRL', 'PEG10-OE-1'), c('OE-CTRL', 'PEG10-OE-1+TGFb1'), c('OE-CTRL', 'PEG10-OE-2'), c('OE-CTRL', 'PEG10-OE-2+TGFb1')), test = "t.test", step_increase = 0.1, textsize = 6, orientation = "x", y_position = 77) + 
  facet_wrap(~variable, nrow = 1, scales = "free_y") + ylab("Number") -> p6

pdf("./Figure9c_9d.pdf", width = 6, height = 16)
(p5/p6)
dev.off() 
