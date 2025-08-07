ml R/4.3.3

# Merging quantification files
merge=$(echo julia src/02_analysis/merge_quantification.jl --dir [path_to_STAR_folder])
merge=$(echo $merge --counts 3 --header 4 --skip_pattern "ENSG" --output STAR_count.csv)
merge=$(echo $merge --gencode 40)
$merge

# DEG analysis
cmd=$(echo Rscript src/02_analysis/deseq.R --counts STAR_counts.csv --metadata metadata.csv)
cmd=$(echo $cmd --comparisons comparisons.txt --sample_col altFilename)
cmd=$(echo $cmd --logFC 1 --use_padj TRUE --padj 0.05 --output_dir .)
$cmd
