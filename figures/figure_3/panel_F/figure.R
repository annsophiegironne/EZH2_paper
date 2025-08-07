library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# To see how DEG analysis was computed, see figures/figure_3/panel_F/commands.sh
data <- read.table("figures/figure_3/panel_F/UNC1999_siATF4_vs_UNC1999_siControl_sorted_logFC_geneID.txt", header=TRUE)

data$deg <- ifelse(data$log2FoldChange>=1&data$pvalue<0.05, "Upregulated",
            ifelse(data$log2FoldChange<=(-1)&data$pvalue<0.05, "Downregulated",
            "Not significant"))

data$labels <- ifelse(data$deg!="Not significant", data$gene_id, NA)
data$labels <- gsub("^ENSG", NA, data$labels)

filename <- "figures/figure_3/panel_F/figure_3_F.pdf"
plot <- ggplot(data, aes(x=log2FoldChange, y=-log10(pvalue), color=deg, label=labels))+
            geom_point(size=1)+
            scale_color_manual(values=c("Not significant"="black",
                                        "Upregulated"="#e69300",
                                        "Downregulated"="#5789c2"))+
            labs(color="Expression")+
            theme_minimal()+
            theme(axis.title=element_text(size=15),
                legend.text=element_text(size=13),
                axis.text=element_text(size=15),
                legend.title=element_text(size=13))+
            geom_hline(yintercept=-log10(0.05), linetype=2)+
            geom_vline(xintercept=c(-1, 1), linetype=2)+
            geom_text_repel(max.overlaps=50, show.legend=FALSE)
ggsave(filename, plot, height=5, width=7)
