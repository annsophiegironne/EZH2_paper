library(fgsea)
library(tidyr)
library(dplyr)
library(ggplot2)
library(argparser)
library(stringr)
library(data.table)

initialize_parser <- function() {
    parser <- arg_parser("Calculate pathway enrichment using fgsea")
    parser <- add_argument(parser, "--file", type="character",
                        help="Path to DESeq2 output file. This file should be unfiltered and still contain NA values.")
    parser <- add_argument(parser, "--pathway", type="character", nargs="*", default="Hallmark",
                        help="Pathway to use for GSEA analysis. One of GOBP, Hallmark or Reactome.")
    parser <- add_argument(parser, "--analysis_name", type="character", default="results",
                        help="Comparison identification, e.g 'treatment_vs_control'")
    parser <- add_argument(parser, "--col", type="character", default="log2FoldChange",
                        help="Column name of log2FoldChange, e.g. log2FoldChange when using Limma")
    parser <- add_argument(parser, "--plot", default=TRUE, type="logical",
                        help="If TRUE, a bar plot will be generated to visualize the top 20 enriched pathways.")
    parser <- add_argument(parser, "--output_dir", default="./", type="character",
                        help="Path to the directory where results will be saved. If the directory does not exist, it will be created")
    return(parser)
}

verifications <- function(args) {
    if (is.na(args$file)) {
        stop("No file given.")
    }

    if (!file.exists(args$file)) {
        stop("Path to file does not exist.")
    }
}

prepare_output_dir <- function(output_dir) {
    if (!file.exists(output_dir)) {
        split <- str_split(output_dir, "/")[[1]]
        output_dir <- split[1]
        for (sub in split[2:length(split)]) {
            output_dir <- paste(output_dir, sub, sep="/")
            if (!file.exists(output_dir)) {
                dir.create(output_dir)
            }
        }
    }

    if (!endsWith(output_dir, "/")) {
        output_dir <- paste0(output_dir, "/")
    }

    return(output_dir)
}


get_pathway <- function(pathway) {
    if (pathway == "GOBP") {
        pathway <- gmtPathways("annotations/GSEA/c5.go.bp.v2023.1.Hs.symbols.gmt")
    } else if (pathway == "Hallmark") {
        pathway <- gmtPathways("annotations/GSEA/h.all.v2023.1.Hs.symbols.gmt")
    } else if (pathway == "Reactome") {
        pathway <- gmtPathways("annotations/GSEA/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
    } else {
        stop("Unrecognized pathway parameter.")
    }
    return(pathway)
}

run_gsea <- function(gene_list, pathway, analysis_name, output_dir, plot) {
    fgseaRes <- fgsea(pathways = pathway,
                stats    = gene_list,
                eps      = 0.0,
                minSize  = 15,
                maxSize  = 500,
                nPermSimple = 100000)
    pathway_name <- tolower(gsub("_.*", "", names(pathway[1])))

    # All output
    filename <- paste0(output_dir, pathway_name, "_", analysis_name, ".txt")
    fwrite(fgseaRes[order(-NES, pval)], filename, sep="\t", sep2=c("", " ", ""), row.names=FALSE)

    # Significant pathways
    filename <- gsub(".txt", "_sig.txt", filename)
    fwrite(fgseaRes[pval<0.05][order(-NES, pval)], filename, sep="\t", sep2=c("", " ", ""), row.names=FALSE)


    if (plot == TRUE) {
        filename <- gsub("sig.txt", "topPathways.pdf", filename)
        up <- head(fgseaRes[pval<0.05][order(abs(NES), decreasing=TRUE), ][NES>0], n=10)
        down <- head(fgseaRes[pval<0.05][order(abs(NES), decreasing=TRUE), ][NES<0], n=10)
        top <- rbind(up, down)

        if (pathway_name == "gobp") {
            width <- 15
        } else {
            width <- 10
        }

        plot <- ggplot(top, aes(y=reorder(pathway, NES), x=NES))+
                        geom_bar(stat="identity", fill="black")+
                        theme_void()+
                        labs(title=pathway_name, y="Pathways")+
                        theme(panel.grid.minor=element_blank(),
                            panel.grid.major=element_line(color="#F2F2F2"),
                            axis.title=element_text(size=17),
                            axis.text=element_text(size=15),
                            axis.title.y=element_blank(),
                            plot.title=element_text(size=17))
        ggsave(filename, plot, width=width, height=7)
    }
}

main <- function(args) {
    output_dir <- prepare_output_dir(args$output_dir)
    pathways <- str_split_1(args$pathway, ",")


    data <- read.table(args$file, header=TRUE)

    if (args$col != "log2FoldChange") {
        names(data)[names(data)==args$col] <- "log2FoldChange"
    }

    data <- data %>% select(log2FoldChange, gene_id) %>% na.omit() %>%
            arrange(log2FoldChange) %>% group_by(log2FoldChange) %>%
            mutate(random_order = runif(n())) %>% arrange(-log2FoldChange, random_order)

    gene_list <- data$log2FoldChange
    names(gene_list) <- data$gene_id

    rank <- data %>% select(gene_id, log2FoldChange)
    write.table(rank, file=file.path(output_dir, paste0(args$analysis_name, "_for_GSEAPreranked.rnk")), 
                col.names=FALSE, quote=FALSE, row.names=FALSE, sep="\t")

    for (pathway in pathways) {
        pathway <- get_pathway(pathway)
        run_gsea(gene_list, pathway, args$analysis_name, output_dir, args$plot)
    }   
}


parser <- initialize_parser()
args <- parse_args(parser)
verifications(args)
main(args)
