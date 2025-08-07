library(dplyr)
library(tidyr)
library(DiffBind)

df <- read.csv("metadata_PDX1915_diffbind.csv", header=TRUE)
counts <- dba(sampleSheet=df, bRemoveRandom=TRUE, minOverlap=1,
            config=data.frame(fragmentSize=0, doBlacklist=FALSE, doGreylist=FALSE))

# See figure_1/panel_F/commands.sh to see how all_peaks_called.bed was generated
H3K27me3.peaks <- read.table("PDX1915_H3K27me3/all_peaks_called.bed", header=FALSE)
counts <- dba.count(counts, peaks=H3K27me3.peaks, bParallel=FALSE, 
                    summits=FALSE, minOverlap=1, bSubControl=TRUE, bScaleControl=TRUE)

dba.save(counts, file="H3K27me3_peaks", pre="dba_", dir="figures/figure_1/panel_H")