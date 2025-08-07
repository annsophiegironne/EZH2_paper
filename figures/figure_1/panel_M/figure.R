library(dplyr)
library(tidyr)
library(ggplot2)

# To see how enrichment scores were computed, see figures/figure_1/panel_L/commands.sh
data <- read.csv("figures/figure_1/panel_L/gene_id_up_enrichmentScores_logFC-1_padj5%.csv", header=T)
data <- data[data$class_id!="DNA"&data$class_id!="RC", ]
data <- data[data$pvalue<0.05, ]

data$class_id <- gsub("Retroposon", "SVA", data$class_id)
class_palette <- c("SVA"="#8da0cb", "LTR"="#fc8d62", "LINE"="#66c2a5", "SINE"="#e78ac3")

filename <- "figures/figure_1/panel_M/figure_1_M.pdf"
plot <- ggplot(data, aes(x=log2OR, fill=class_id))+
            geom_density()+
            theme_bw()+
            theme(panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="#F2F2F2"),
                axis.title=element_text(size=15),
                axis.text=element_text(size=14),
                legend.title=element_text(size=14),
                legend.text=element_text(size=13),
                strip.background=element_blank(),
                strip.text.y=element_text(size=15))+
            geom_vline(xintercept=c(-1,1), linetype="dashed", color="black")+
            labs(fill="TE class", y="Density")+
            scale_fill_manual(values=class_palette)+
            facet_grid(rows=vars(class_id), scales="free")+
            xlim(-5, 5)
ggsave(filename, plot)
