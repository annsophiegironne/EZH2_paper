library(dplyr)
library(tidyr)
library(ggplot2)

# To see how DEG analysis was computed, see figures/figure_3/panel_F/commands.sh
unc <- read.table("figures/figure_3/panel_H/UNC1999_siControl_vs_DMSO_siControl_sorted_logFC_geneID.txt", header=TRUE)
unc <- unc %>% select(-baseMean, -lfcSE, -stat, -padj)
names(unc) <- c("gene_id", "logFC_UNC1999", "pvalue_UNC1999")

si <- read.table("UNC1999_siATF4_vs_UNC1999_siControl_sorted_logFC_geneID.txt", header=TRUE)
si <- si %>% select(-baseMean, -lfcSE, -stat, -padj)
names(si) <- c("gene_id", "logFC_siATF4", "pvalue_siATF4")

data <- inner_join(unc, si, by="gene_id")
data$color <- ifelse(data$logFC_siATF4>1&abs(data$logFC_UNC1999)>1, "Upregulated",
              ifelse(data$logFC_siATF4<(-1)&abs(data$logFC_UNC1999)>1, "Downregulated",
              "Not significant"))

pearson <- cor.test(data$logFC_siATF4, data$logFC_UNC1999, method="pearson")
corr <- pearson$estimate[[1]]
pvalue <- pearson$p.value

if (pvalue<0.0001) {
    pvalue.text <- "P < 0.0001"
} else {
    pvalue.text <- paste0("P = ", round(pvalue,3))
}

filename <- "figures/figure_3/panel_H/figure_3_H.pdf"
plot <- ggplot(data, aes(x=logFC_siATF4, y=logFC_UNC1999, color=color))+
            geom_point(size=0.5)+
            theme_bw()+
            labs(x="Log2FoldChange\n(siATF4 vs siControl in UNC1999)",
                 y="Log2FoldChange\n(UNC1999 vs DMSO with siControl)",
                 color="Expression")+
            scale_color_manual(values=c("Upregulated"="#e69300",
                                        "Downregulated"="#5789c2",
                                        "Not significant"="#919294"))+
            theme(axis.text=element_text(size=15),
                  legend.title=element_text(size=14),
                  legend.text=element_text(size=13))+
            annotate("text", x=3, y=4, size=3, label=paste0("R = ", round(corr, 3)))+
            annotate("text", x=3, y=3.7, size=3, label=pvalue.text)+
            geom_hline(yintercept=0, color="black", linewidth=0.01)+
            geom_vline(xintercept=0, color="black", linewidth=0.01)+
            scale_x_continuous(breaks=seq(-6, 6, 2), limits=c(-6, 6))+
            scale_y_continuous(breaks=seq(-6, 6, 1), limits=c(-6, 6))
ggsave(filename, plot, height=4, width=7)
