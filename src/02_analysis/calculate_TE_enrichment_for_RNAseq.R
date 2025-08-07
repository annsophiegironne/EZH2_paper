library(argparser)
library(dplyr)
library(stringr)

initialize_parser <- function() {
    parser <- arg_parser("Calculate TE enrichment for subfamilies/families/classes in given list of TEs")
    parser <- add_argument(parser, "--file", type="character", help="Path to the unfiltered DESeq2 results for TEs. 
                                                                    This file must contain log2FoldChange, pvalue and/or padj columns, as well as the gene_id column with TE ids.")
    parser <- add_argument(parser, "--logFC", default=0.3, type="numeric",
                            help="Threshold for log2 fold change (logFC) to identify differentially expressed TEs")
    parser <- add_argument(parser, "--use_padj", default=TRUE, type="logical",
                            help="If TRUE, the adjusted p-value (padj) will be used to determine statistical significance. If FALSE, the pvalue will be used instead")
    parser <- add_argument(parser, "--padj", default=0.05, type="numeric",
                            help="Threshold for adjusted p-value (padj) (or pvalue) to determine statistical significance")
    parser <- add_argument(parser, "--in_level", default="copy_id", type="character",
                            help="Specify the level of TE classification on which DEG analysis was computed:\n
                                - 'copy_id' : DEG analysis was computed on TE insertions\n
                                - 'gene_id' : DEG analysis was computed on TE subfamilies (TE insertions of each subfamily were summed)")
    parser <- add_argument(parser, "--level", default="gene_id", type="character",
                            help="Specify the level of analysis for the TE enrichment:\n
                                - 'gene_id': analyse at the subfamily level (ie. enrichment scores are calculated based on differentially expressed (in_level) for each subfamily)\n
                                - 'family_id': analyse at the family level (ie. enrichment scores are calculated based on differentially expressed (in_level) for each family)\n
                                - 'class_id' : analyse at the class level (ie. enrichment scores are calcualted based on differentially expressed (in_level) for each class)")
    parser <- add_argument(parser, "--output_dir", default="./", type="character",
                            help="Path to the directory where results will be saved. If the directory does not exist, it will be created")
    return(parser)
}


verify_args <- function(args) {
    if (is.na(args$file)) {
        stop("No file given")
    }

    if (!(args$in_level %in% c("copy_id", "gene_id"))) {
        stop("--in_level not recognized")
    }

    if (!(args$level %in% c("gene_id", "family_id", "class_id"))) {
        stop("--level not recognized")
    }

    if (!file.exists(args$output_dir)) {
        dir.create(args$output_dir)
    }
}


compute_fisher <- function(a, b, c, d) {
    matrix <- matrix(c(a, c, b, d), nrow=2, ncol=2)
    results <- fisher.test(matrix, conf.int=TRUE, conf.level=0.95, alternative="two.sided")
    return(results)
}

extract_results <- function(results) {
    data <- c(log2(results$estimate[[1]]), results$p.value, 
            log2(results$conf.int[1]), log2(results$conf.int[2]))
    return(data)
}

read_annotations <- function() {
    gtf <- read.csv("annotations/TE/TE_Homo_sapiens.GRCh38.annotations.csv", header=T, sep=",")
    return(gtf)
}

create_dataframe <- function(level, gtf) {
    size <- length(unique(gtf[[level]]))
    df <- setNames(as.data.frame(matrix(nrow=size, ncol=6)),
                c(level, "class_id", "log2OR", "pvalue", "conf_lower", "conf_upper"))

    return(df)
}

iterate_and_store_results <- function(df, file, in_level, level, gtf) {
    i <- 1
    for (data in unique(gtf[[level]])) {
        abcd <- compute_cells(data, in_level, level, file, gtf)
        results <- compute_fisher(abcd[1], abcd[2], abcd[3], abcd[4])
        results <- extract_results(results)
        results <- c(data, gtf[gtf[[level]]==data,]$class_id[1], results)
        df[i, ] <- results
        i <- i + 1
    }
    return(df)
}

compute_cells <- function(data, in_level, level, file, gtf) {
    # In level and in DEG
    a <- sum(file[[level]] == data)

    # In level but not in DEG
    b <- sum(gtf[[level]] == data) - a

    # In DEG but not in level
    c <- sum(file[[level]] != data)

    # Not in DEG and not in level
    d <- length(gtf[[in_level]]) - a - b - c
    return(c(a, b, c, d))
}

filter_gtf <- function(gtf, in_level, level) {
    return(unique(gtf[,c(in_level, level, "class_id")]))
}


write_results <- function(df, level, deg, logFC, use_padj, padj, output_dir) {
    if (!endsWith(output_dir, "/")) {
        output_dir <- paste0(output_dir, "/")
    }
    if (use_padj == TRUE) {
        pval <- "padj"
    } else {
        pval <- "pvalue"
    }

    treshold <- padj * 100

    filename <- paste0(output_dir, level, "_", deg, "_enrichmentScores_logFC-", logFC, "_", pval, treshold, "%.csv")
    write.csv(df, file=filename, quote=FALSE, row.names=FALSE)
}

filter_file <- function(file, logFC, use_padj, padj, in_level) {
    if (use_padj == TRUE) {
        col <- "padj" 
    } else {
        col <- "pvalue"
    }

    file <- na.omit(subset(file, select=c(in_level, "log2FoldChange", col))[file[[col]]<=padj,])
    files <- list()

    files$up <- file %>% filter(log2FoldChange >= abs(logFC))
    files$down <- file %>% filter(log2FoldChange <= -(abs(logFC)))

    return(files)
}


change_ids <- function(file, in_level) {
    if (in_level == "copy_id") {
        file$gene_id <- gsub(":.*", "", file$gene_id)
    } else {
        file$gene_id <- gsub("(_dup.*)*:.*", "", file$gene_id)
    }

    colnames(file)[colnames(file)=="gene_id"] <- in_level

    return(file)
}

main <- function(file, logFC, use_padj, padj, in_level, level, output_dir) {
    file <- read.table(file, header=TRUE, sep="\t")

    gtf <- read_annotations()
    gtf <- filter_gtf(gtf, in_level, level)

    file <- change_ids(file, in_level)
    files <- filter_file(file, logFC, use_padj, padj, in_level)

    for (deg in names(files)) {
        data <- inner_join(files[[deg]], gtf, by=in_level)
        df <- create_dataframe(level, gtf)
        df <- iterate_and_store_results(df, data, in_level, level, gtf)
        write_results(df, level, deg, logFC, use_padj, padj, output_dir)
    }

}

parser <- initialize_parser()
args <- parse_args(parser)
verify_args(args)
main(args$file, args$logFC, args$use_padj, args$padj, args$in_level, args$level, args$output_dir)