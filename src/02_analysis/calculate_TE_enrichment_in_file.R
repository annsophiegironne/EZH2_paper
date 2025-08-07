library(argparser)
library(dplyr)

initialize_parser <- function() {
    parser <- arg_parser("Calculate TE enrichment for subfamilies/families/classes in given list of TEs")
    parser <- add_argument(parser, "--file",
                        help="BED/CSV file of TEs")
    parser <- add_argument(parser, "--header", default=FALSE,
                        help="Whether file contains a header row")
    parser <- add_argument(parser, "--level", default="subfamilies",
                        help="One of gene_id (subfamily), family_id (family) or class_id (class)")
    parser <- add_argument(parser, "--in_level", default="insertions",
                        help="One of copy_id (insertion), gene_id (subfamily)")
    parser <- add_argument(parser, "--col", default=1,
                        help="Column in file where TE ids are")
    parser <- add_argument(parser, "--output_dir", default="./",
                        help="Where to store the output file")
    return(parser)
}

verify_args <- function(args) {
    if (!(args$in_level %in% c("copy_id", "gene_id"))) {
        stop("--in_level not recognized")
    }

    if (!(args$level %in% c("gene_id", "family_id", "class_id"))) {
        stop("--level not recognized")
    }

    if (args$col < 1) {
        stop("Invalid column number")
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
            log2(results$conf.int[1]), log2(results$conf.int[2]+0.001))
    return(data)
}

read_annotations <- function() {
    gtf <- read.csv("annotations/TE/TE_Homo_sapiens.GRCh38.annotations.csv", header=TRUE, sep=",")
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
        results <- c(data, gtf[gtf[[level]]==data, ]$class_id[1], results)
        df[i, ] <- results
        i <- i + 1
    }
    return(df)
}

compute_cells <- function(data, in_level, level, file, gtf) {
    a <- sum(file[[level]] == data)
    b <- sum(gtf[[level]] == data) - a
    c <- sum(file[[level]] != data)
    d <- length(gtf[[in_level]]) - a - b - c
    return(c(a, b, c, d))
}


filter_gtf <- function(gtf, in_level, level) {
    return(unique(gtf[,c(in_level, level, "class_id")]))
}

write_results <- function(df, level, output_dir) {
    if (!endsWith(output_dir, "/")) {
        output_dir <- paste0(output_dir, "/")
    }
    filename <- paste0(output_dir, level, "_enrichmentScores.csv")
    write.csv(df, file=filename, quote=FALSE, row.names=FALSE)
}

change_ids <- function(file, in_level, col) {
    colnames(file)[col] <- in_level

    if (in_level == "copy_id") {
        file$copy_id <- gsub(":.*", "", file$copy_id)
    } else if (in_level == "gene_id") {
        file$gene_id <- gsub("(_dup.*)*:.*", "", file$gene_id)
    }

    return(file)
}

main <- function(args) {
    file <- read.table(args$file, header=args$header, sep="\t")
    
    gtf <- read_annotations()
    gtf <- filter_gtf(gtf, args$in_level, args$level)


    file <- change_ids(file, args$in_level, args$col)
    file <- inner_join(file, gtf, by=args$in_level)

    df <- create_dataframe(args$level, gtf)
    df <- iterate_and_store_results(df, file, args$in_level, args$level, gtf)
    write_results(df, args$level, args$output_dir)
}

parser <- initialize_parser()
args <- parse_args(parser)
verify_args(args)
main(args)