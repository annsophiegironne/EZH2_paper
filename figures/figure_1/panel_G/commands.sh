ml deeptools/3.5.1

# 'all_peaks.bed' are all peaks (common, CTL, UNC1999) concatenated after DiffBind analysis
# Each bdg sample file was scaled, then the mean was computed on all bdg scaled files
# then BedGraphToBigWig (ucsc) was used to convert to bigwig file
# eg:
# julia src/02_analysis/scale_bdg.jl --file file1.bw --factor 1.2304 --output file1.scaled.bdg
# julia src/02_analysis/compute_mean_bdg.jl --files file1.bw file2.bw file3.bw --output files.mean.bdg

# Heatmap
matrix=$(echo computeMatrix reference-point --referencePoint center -a 2000 -b 2000)
matrix=$(echo $matrix -R specific_to_CTL.bed common.bed specific_to_UNC1999.bed)
matrix=$(echo $matrix -S CTL.mean.bw UNC1999.mean.bw --skipZeros)
matrix=$(echo $matrix --missingDataAsZero -p max -o diff_peaks.mat.gz)
$matrix

heatmap=$(echo plotHeatmap -m diff_peaks.mat.gz -o diff_peaks.mat-heatmap.pdf --dpi 500)
heatmap=$(echo $heatmap --heatmapHeight 3 --heatmapWidth 3 --colorMap Oranges)
$heatmap

# Profile
matrix=$(echo computeMatrix reference-point --referencePoint center -a 2000 -b 2000)
matrix=$(echo $matrix -R all_peaks.bed)
matrix=$(echo $matrix -S CTL.mean.bw UNC1999.mean.bw --skipZeros)
matrix=$(echo $matrix --missingDataAsZero -p max -o all_peaks.mat.gz)
$matrix

heatmap=$(echo plotHeatmap -m all_peaks.mat.gz -o all_peaks.mat-heatmap.pdf --dpi 500)
heatmap=$(echo $heatmap --heatmapHeight 3 --heatmapWidth 3 --colorMap Oranges)
$heatmap

profile=$(echo plotProfile -m all_peaks.mat.gz -o all_peaks.mat-heatmap.pdf --dpi 500)
$profile


