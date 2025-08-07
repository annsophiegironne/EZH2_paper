library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

k27ac.dmso <- read.delim("figures/figure_S1/panel_M/common_and_specific_to_DMSO_annotations.txt", header=T, sep="\t")
k27ac.unc <- read.delim("figures/figure_S1/panel_M/specific_to_UNC1999_annotations.txt", header=T, sep="\t")
k27me3 <- read.delim("figures/figure_S1/panel_M/all_peaks_called_annotations.txt", header=T, sep="\t")

k27ac.dmso$condition <- "DMSO"
k27ac.unc$condition <- "UNC"
k27me3$condition <- "H3K27me3"

k27 <- rbind(k27ac.dmso, k27ac.unc, k27me3)
k27 <- k27 %>% select(condition, annotation)

k27$myanno <- ifelse(grepl("Promoter", k27$annotation), "Promoter",
              ifelse(grepl("Intergenic|Downstream", k27$annotation), "Distal intergenic",
              ifelse(grepl("Exon|UTR", k27$annotation), "Exon/UTR",
              ifelse(grepl("Intron", k27$annotation), "Intron", ""))))

# Verifications
dim(k27[k27$myanno=="",])[1] == 0
dim(k27[k27$myanno=="Promoter",])[1] == dim(k27[grepl("Promoter", k27$annotation),])[1]
dim(k27[k27$myanno=="Distal intergenic",])[1] == dim(k27[grepl("Intergenic|Downstream", k27$annotation),])[1]
dim(k27[k27$myanno=="Intron",])[1] == dim(k27[grepl("Intron", k27$annotation),])[1]
dim(k27[k27$myanno=="Exon/UTR",])[1] == dim(k27[grepl("Exon|UTR", k27$annotation),])[1]


data <- k27 %>% group_by(condition, myanno) %>% summarize(n=n()) %>% mutate(perc = round(n/sum(n)*100,1))
data$condition <- factor(data$condition, levels=c("H3K27me3", "UNC", "DMSO"))
data$myanno <- factor(data$myanno, levels=c("Distal intergenic", "Intron", "Exon/UTR", "Promoter"))

filename <- "figures/figure_S1/panel_M/figure_S1_M.pdf"
plot <- ggplot(data, aes(x=perc, y=condition, fill=myanno))+
            geom_col()+
            theme_bw()+
            theme(axis.title=element_text(size=15),
                axis.text=element_text(size=14),
                legend.title=element_text(size=14),
                legend.text=element_text(size=14))+
            labs(fill="Annotations", x="% loci distribution", y="Condition")+
            scale_fill_brewer(palette="Paired")
ggsave(filename, plot, width=10, height=5)
