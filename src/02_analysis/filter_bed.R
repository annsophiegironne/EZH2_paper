library(argparser)
library(stringr)

initialize_parser <- function() {
    parser <- arg_parser("Filter differential binding analysis from DiffBind")
    parser <- add_argument(parser, "--file", type="character",
                        help="Path to binding analysis file. This file should include all peaks as well as data from the differential binding analysis.
                        Usually a {analysis}_wholeReport.txt file")
    parser <- add_argument(parser, "--comparison", type="character",
                        help="Names of the two groups to compare (e.g., 'treatment,control'). Peaks are classified as specific to treatment
                        (logFC > treshold and pvalue < treshold), specific to control (logFC < -treshold, pvalue < treshold) or common (all others)")
    parser <- add_argument(parser, "--logFC", default=1, type="numeric",
                        help="Threshold for log2 fold change (logFC) to identify differentially bound peaks")
    parser <- add_argument(parser, "--use_padj", default=FALSE, type="logical",
                        help="If TRUE, the adjusted p-value (padj) will be used to determine statistical significance. If FALSE, the pvalue will be used instead")
    parser <- add_argument(parser, "--padj", default=0.05, type="numeric",
                        help="Threshold for adjusted p-value (padj) (or pvalue) to determine statistical significance")
    parser <- add_argument(parser, "--output_dir", default=".", type="character",
                        help="Path to the directory where results will be saved. If the directory does not exist, it will be created")
    return(parser)
}

verifications <- function(args) {
    if (is.na(args$file)) {
        stop("No file given.")
    }

    if (args$padj < 0 || args$padj > 1) {
        stop("padj/pvalue must be between 0 and 1.")
    }

    if (is.na(args$comparison)) {
        stop("No comparison given.")
    }
}

create_output_dir <- function(output_dir) {
    if (!file.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }

}


read_file <- function(file) {
    split <- str_split(file, fixed("."))[[1]]
    ext <- split[length(split)]

    if (ext == "csv") {
        file <- read.csv(file, header=TRUE, check.names=FALSE)
    } else if (ext == "txt" || ext == "tsv") {
        file <- read.table(file, header=TRUE, sep="\t", check.names=FALSE)
    } else {
        stop("File must be in .csv, .txt or .tsv format.")
    }
    return(file)
}


identify_DBPs <- function(data, treatment, control, logFC, use_padj, padj, output_dir) {
    if (use_padj == TRUE) {
        pval <- "FDR"
    } else {
        pval <- "p-value"
    }

    trmt <- data[data[[pval]]<=padj&data$Fold>=logFC, ]
    filename <- file.path(output_dir, paste0("specific_to_", treatment, ".txt"))
    write.table(trmt, file=filename, col.names=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
    filename <- gsub(".txt", ".bed", filename)
    write.table(trmt[,1:3], file=filename, col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")

    ctrl <- data[data[[pval]]<=padj&data$Fold<=(-logFC), ]
    filename <- file.path(output_dir, paste0("specific_to_", control, ".txt"))
    write.table(ctrl, file=filename, col.names=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
    filename <- gsub(".txt", ".bed", filename)
    write.table(ctrl[,1:3], file=filename, col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")

    common <- data[(abs(data$Fold)<=logFC)|(data[[pval]]>padj&data$Fold>=abs(logFC)), ]
    filename <- file.path(output_dir, paste0("common.txt"))
    write.table(common, file=filename, col.names=TRUE, quote=FALSE, row.names=FALSE, sep="\t")
    filename <- gsub(".txt", ".bed", filename)
    write.table(common[,1:3], file=filename, col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")
}


main <- function(args) {
    create_output_dir(args$output_dir)
    file <- read_file(args$file)
    comp <- str_split_1(args$comparison, ",")
    treatment <- comp[1]
    control <- comp[2]

    identify_DBPs(file, treatment, control, args$logFC, args$use_padj, args$padj, args$output_dir)
}


parser <- initialize_parser()
args <- parse_args(parser)
verifications(args)
main(args)