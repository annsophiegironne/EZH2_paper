ml deeptools/3.5.1

# 'all_peaks_called.bed' are peaks called in CTL/UNC1999 samples from MACS3 concatenated
# and merged with bedtools merge (v2.30.0)
# e.g.:
# cat MDAMB436_H3K27me3*.bed > all_peaks.bed
# sort -k1,1 -k2,2n all_peaks.bed > all_peaks.bed.tmp
# bedtools merge -i all_peaks.bed.tmp > all_peaks.bed
#
# Each bdg sample file was scaled, then the mean was computed on all bdg scaled files
# then BedGraphToBigWig (ucsc) was used to convert to bigwig file
# eg:
# julia src/02_analysis/scale_bdg.jl --file file1.bw --factor 1.2304 --output file1.scaled.bdg
# julia src/02_analysis/compute_mean_bdg.jl --files file1.bw file2.bw file3.bw --output files.mean.bdg


matrix=$(echo computeMatrix reference-point --referencePoint center -a 2000 -b 2000)
matrix=$(echo $matrix -R all_peaks_called.bed -S CTL.mean.bw UNC1999.mean.bw --skipZeros)
matrix=$(echo $matrix --missingDataAsZero -p max -o all_peaks_called.mat.gz)
$matrix

heatmap=$(echo plotHeatmap -m all_peaks_called.mat.gz -o all_peaks_called.mat-heatmap.pdf --dpi 500)
heatmap=$(echo $heatmap --heatmapHeight 3 --heatmapWidth 3 --colorMap Blues)
$heatmap