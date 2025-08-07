ml bedtools/2.30.0

# left_upper_corner.bed was generated from figures/figure_1/panel_H/figure.R with
# corner <- data[data$H3K27ac>0&data$H3K27me3<0,]

# Overlap
bedtools intersect -a left_upper_corner_H3K27me3.bed -b TE.GRCh38.bed -F 1 -wb | cut -f 4,5,6,7 > TEs_fully_in_left_upper_corner.bed

# Enrichment values
cmd=$(echo Rscript src/02_analysis/calculate_TE_enrichment_in_file.R)
cmd=$(echo $cmd --file TEs_fully_in_left_upper_corner.bed --col 4)
$cmd