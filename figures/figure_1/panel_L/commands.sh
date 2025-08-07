ml R/4.3.3

# Merging quantification files
merge=$(echo julia src/02_analysis/merge_quantification.jl --dir [path_to_salmon_folder] --pattern quant.sf)
merge=$(echo $merge --counts 5 --header 1 --skip_pattern "ENST|ENSMUST|_MM" --output salmon_count.csv)
merge=$(echo $merge --gencode 40)
$merge

# DEG analysis
cmd=$(echo Rscript src/02_analysis/deseq.R --counts salmon_counts.csv --metadata metadata.csv)
cmd=$(echo $cmd --comparisons comparisons.txt --sample_col altFilename --TE TRUE --level insertion)
cmd=$(echo $cmd --logFC 1 --use_padj TRUE --padj 0.05 --output_dir .)
$cmd

# Enrichment values
cmd=$(echo Rscript src/02_analysis/calculate_TE_enrichment_for_RNAseq.R)
cmd=$(echo $cmd --file UNC1999_vs_CTL_unfiltered.txt --logFC 1 --use_padj TRUE)
cmd=$(echo $cmd --padj 0.05 --level gene_id --in_level copy_id --output_dir figures/figure_1/panel_L)
$cmd