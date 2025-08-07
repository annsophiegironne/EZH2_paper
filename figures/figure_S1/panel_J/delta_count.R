library(dplyr)
library(tidyr)
library(DiffBind)

df <- read.csv("metadata_PDX1915_diffbind.csv", header=TRUE)
counts <- dba(sampleSheet=df, bRemoveRandom=TRUE, minOverlap=1,
            config=data.frame(fragmentSize=0, doBlacklist=FALSE, doGreylist=FALSE))

H3K27ac.peaks <- read.table("PDX1915_H3K27ac/specific_to_UNC1999.bed", header=FALSE)
counts <- dba.count(counts, peaks=H3K27ac.peaks, bParallel=FALSE,
                    summits=FALSE, minOverlap=1, bSubControl=TRUE, bScaleControl=TRUE)

dba.save(counts, file="H3K27ac_peaks", pre="dba_", dir="figures/figure_S1/panel_J")
