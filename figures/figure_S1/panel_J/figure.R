library(dplyr)
library(tidyr)
library(ggplot2)
library(DiffBind)

# See figures/figure_S1/panel_J/delta_count.R to see how counts were generated
counts <- dba.load(file="H3K27ac_peaks", dir="figures/figure_S1/panel_J")

# See figures/figure_S1/panel_J/commands.R to see how normFactors were generated
H3K27ac.factors <- read.table("figures/figure_S1/panel_J/H3K27ac_normFactors.txt", header=TRUE)
H3K27me3.factors <- read.table("figures/figure_S1/panel_J/H3K27me3_normFactors.txt", header=TRUE)
order <- dba.show(counts)$ID

norm <- rbind(H3K27ac.factors, H3K27me3.factors)
row.names(norm) <- norm$ID
norm.factors <- as.numeric(norm[order, ]$norm.factors)

counts <- dba.normalize(counts, method=DBA_DESEQ2, normalize=norm.factors)
counts <- dba.count(counts, peaks=NULL, score=DBA_SCORE_NORMALIZED)

data <- dba.peakset(counts, bRetrieve=TRUE, DataType=DBA_DATA_FRAME)

data$chr <- paste0(data$CHR, ":", data$START, "-", data$END)
data <- data[,-c(1:3)]
data <- data %>% pivot_longer(cols=!chr, names_to="samples", values_to="counts")
data$condition <- gsub("PDX1915_(.*)_H3K27.*", "\\1", data$samples)
data$mark <- gsub(".*_(H3K27.*)_.*", "\\1", data$samples)
data <- data %>% select(-samples) %>%
        group_by(chr, mark, condition) %>%
        summarize_all(mean)
data <- data %>% spread(condition, counts)
data <- data %>% group_by(chr, mark) %>% mutate(delta = UNC1999 - CTL)
data <- data %>% select(-CTL, -UNC1999) %>% spread(mark, delta)
data <- data[abs(data$H3K27ac)<=600&abs(data$H3K27me3)<=200, ]
data <- data[data$H3K27ac!=0&data$H3K27me3!=0, ]

result <- cor.test(data$H3K27ac, data$H3K27me3, method="pearson")

if (result$p.value < 0.001) {
    pvalue <- "P < 0.001"
} else {
    pvalue <- paste0("P = ", round(result$pvalue, 3))
}

corr <- round(result$estimate[[1]], 4)

filename <- "figures/figure_S1/panel_J/figure_S1_J.pdf"
plot <- ggplot(data, aes(x=H3K27me3, y=H3K27ac))+
                geom_point(size=0.1)+
                theme_bw()+
                ylim(-300, 600)+
                xlim(-200, 200)+
                labs(x=expression(Delta * " H3K27me3"), 
                y=expression(Delta * " H3K27ac"))+
                theme(panel.grid.minor=element_blank(), 
                      panel.grid.major=element_line(color="#F2F2F2"),
                      axis.title=element_text(size=13),
                      axis.text=element_text(size=12))+
                annotate("text", x=100, y=400, label=pvalue)+
                annotate("text", x=100, y=370, label=paste0("R = ", corr))+
                geom_hline(yintercept=0)+
                geom_vline(xintercept=0)
ggsave(filename, plot)
