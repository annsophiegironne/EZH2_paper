library(dplyr)
library(tidyr)

# To see how enrichment scores were computed, see figures/figure_S1/panel_Q/commands.sh
data <- read.csv("figures/figure_S1/panel_Q/gene_id_enrichmentScores.csv", header=TRUE)
gtf <- read.csv("/u/gironnea/my_libraries/annotations/TE/TE_HomoSapiens_GRCh38_annotations.csv", header=TRUE)
gtf <- gtf %>% select(-copy_id) %>% unique()

data <- inner_join(data, gtf, by=c("gene_id", "class_id"))

data$class_id <- ifelse(data$class_id=="Retroposon", "SVA", data$class_id)
data <- data[data$class_id %in% c("LINE", "SINE", "LTR", "SVA"), ]
data$color <- ifelse(abs(data$log2OR)>1&data$pvalue<0.05, data$family_id, "Not significant")

palette <- c("Alu"="#e78ac3", "MIR"="#D352A1",
            "ERV1"="#fc8d62", "ERVK"="#FF692F",
            "ERVL"="#CF7C0E",
            "ERVL-MaLR"="#E8A348",
            "Gypsy"="#FDBF6F",
            "L1"="#66c2a5",
            "L2"="#33A02C",
            "CR1"="#155E0F",
            "RTE-X"="#84D57D",
            "SVA"="#8da0cb")

data$color <- factor(data$color, levels=c("Alu", "MIR", "ERVK", "ERV1", "ERVL", "ERVL-MaLR",
                    "Gypsy", "L1", "L2", "CR1", "RTE-X", "SVA", "Not significant"))

filename <- "figures/figure_S1/panel_Q/figure_S1_Q.pdf"
plot <- ggplot(data, aes(x=log2OR, y=-log10(pvalue), color=color))+
            geom_point(size=2)+
            theme_bw()+
            theme(panel.grid.minor=element_blank(),
                panel.grid.major=element_line(color="#F2F2F2"),
                axis.title=element_text(size=15),
                axis.text=element_text(size=14),
                legend.title=element_text(size=14),
                legend.text=element_text(size=13))+
            geom_hline(yintercept=c(-log10(0.05)), linetype="dashed", color="black")+
            geom_vline(xintercept=c(-1, 1), linetype="dashed", color="black")+
            labs(color="Enrichment")+
            scale_color_manual(values=c(palette, "Not significant"="#6C6C6C"))+
            xlim(-5, 5)
ggsave(filename, plot, width=7, height=4)
