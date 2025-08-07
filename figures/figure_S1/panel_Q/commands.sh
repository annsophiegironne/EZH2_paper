ml bedtools/2.30.0

# Overlap
bedtools intersect -a all_peaks_called.bed -b TE.GRCh38.bed -F 1 -wb | cut -f 4,5,6,7 > TEs_fully_in_H3K27me3_peaks.bed

# Enrichment values
cmd=$(echo Rscript src/02_analysis/calculate_TE_enrichment_in_file.R)
cmd=$(echo $cmd --file TEs_fully_in_H3K27me3_peaks --col 4)
$cmd