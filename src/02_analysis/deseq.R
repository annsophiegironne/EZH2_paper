library(argparser)
library(dplyr)
library(stringr)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(pals)
library(processx)
library(tools)

initialize_parser <- function() {
    parser <- arg_parser("Calculate differential expression analysis using DESeq2")
    parser <- add_argument(parser, "--counts", type="character",
                        help="Path to counts file. This file should include raw count data for all samples that will be used for normalization")
    parser <- add_argument(parser, "--metadata", type="character",
                        help="Path to metadata file. This file should include all the samples that will be used for normalization and it must include a column named 'condition' that contains the experimental conditions (e.g., 'treatment1', 'control1', etc.).
                            The values in the 'condition' column must match the conditions specified in the --comparisons")
    parser <- add_argument(parser, "--sample_col", type="character",
                        help="Name of the column in the metadata file that contains the sample identifiers.
                            These sample names must match exactly those in the counts file to ensure correct mapping.")
    parser <- add_argument(parser, "--condition_col", type="character", default="condition",
                        help="Name of the column in the metadata file that contains the conditions to test.
                            Samples will be compared based on these conditions.")
    parser <- add_argument(parser, "--design", type="character", default="condition",
                        help="Experimental design of the analysis. When adding multiple terms, specify the condition using '',
                        eg. 'condition + model'.")
    parser <- add_argument(parser, "--comparisons", type="character",
                        help="Path to a file containing a list of comparisons to be performed. Each comparison should be on a separate line, with the treatment and control conditions separated by a comma (e.g., 'treatment1,control1').
                            The conditions specified in the file must match the values present in the 'condition' column of the metadata file. Each line represents a comparison to be processed by the script.")
    parser <- add_argument(parser, "--logFC", default=1, type="numeric",
                        help="Threshold for log2 fold change (logFC) to identify differentially expressed genes")
    parser <- add_argument(parser, "--use_padj", default=TRUE, type="logical",
                        help="If TRUE, the adjusted p-value (padj) will be used to determine statistical significance. If FALSE, the pvalue will be used instead")
    parser <- add_argument(parser, "--padj", default=0.05, type="numeric",
                        help="Threshold for adjusted p-value (padj) (or pvalue) to determine statistical significance")
    parser <- add_argument(parser, "--TE", default=FALSE, type="logical",
                        help="If TRUE, the parameter 'level' will specify the TE level (e.g., subfamily, insertion) for analysis.
                            If FALSE, the analysis will focus on gene-level counts instead of TE-level counts.")
    parser <- add_argument(parser, "--level", default="subfamily", type="character",
                        help="Specify the level of analysis for the TE counts:\n
                            - 'subfamily': analyse at the subfamily level (ie. insertions counts of a subfamily are summed)\n
                            - 'insertion': analyse at the insertion level (ie. each individual insertion are analyzed)")
    parser <- add_argument(parser, "--save_counts", default=TRUE, type="logical",
                        help="If TRUE, log2-normalized counts matrix will be saved in the specified output directory")
    parser <- add_argument(parser, "--deg_analysis", default=TRUE, type="logical",
                        help="If TRUE, DEG analysis will be performed. Setting to FALSE allows for plotting purposes without doing a DEG analysis.")
    parser <- add_argument(parser, "--run_GSEA", default=TRUE, type="logical",
                        help="If TRUE, GSEA analysis will be run based on given pathway parameter.")
    parser <- add_argument(parser, "--pathway", type="character", nargs="*", default="Hallmark",
                        help="Pathway to use for GSEA analysis. One of GOBP, Hallmark or Reactome.")
    parser <- add_argument(parser, "--plot_GSEA", type="logical", default=TRUE,
                        help="If TRUE, a bar plot will be generated to visualize the top 20 enriched pathways.")
    parser <- add_argument(parser, "--plot_PCA", default=TRUE, type="logical",
                        help="If TRUE, a PCA plot will be generated to visualize sample clustering")
    parser <- add_argument(parser, "--plot_volcano", default=TRUE, type="logical",
                        help="If TRUE, a volcano plot will be generated to visualize differential gene expression data")
    parser <- add_argument(parser, "--plot_sizeFactors", default=FALSE, type="logical",
                        help="If TRUE, size factors plot will be generated to compare normalization factors used across samples")
    parser <- add_argument(parser, "--plot_clustering", default=TRUE, type="logical",
                        help="If TRUE, clustering heatmap will be generated to compare Pearson correlations between samples")
    parser <- add_argument(parser, "--output_dir", default="./", type="character",
                        help="Path to the directory where results will be saved. If the directory does not exist, it will be created")
    
    return(parser)
}

verifications <- function(args) {
    if (is.na(args$counts)) {
        stop("No counts file given.")
    }

    if (!file.exists(args$counts)) {
        stop("Path to counts file does not exist.")
    }

    if (is.na(args$metadata)) {
        stop("No metadata file given.")
    }

    if (!file.exists(args$metadata)) {
        stop("Path to counts file does not exist.")
    }

    if (is.na(args$comparisons) && args$deg_analysis == TRUE) {
        stop("No comparisons file given.")
    }

    if (!is.na(args$comparisons) && !file.exists(args$comparisons)) {
        stop("Path to comparison file does not exist.")
    }

    if (args$padj < 0 || args$padj > 1) {
        stop("padj must be between 0 and 1.")
    }

    if (!(args$level %in% c("subfamily", "insertion"))) {
        stop("Level parameter must be set to either 'subfamily' or 'insertion'")
    }

    if (is.na(args$sample_col)) {
        stop("Sample_col cannot be empty.")
    }
}

create_output_dir <- function(output_dir) {
    if (!file.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }

}

read_counts <- function(counts) {
    ext <- file_ext(counts)

    if (ext == "csv") {
        counts <- read.csv(counts, header=TRUE, check.names=FALSE)
    } else if (ext == "txt" || ext == "tsv") {
        counts <- read.table(counts, header=TRUE, sep="\t", check.names=FALSE)
    } else {
        stop("Counts file must be in .csv, .txt or .tsv format.")
    }

    return(counts)
}

read_metadata <- function(metadata) {
    ext <- file_ext(metadata)

    if (ext == "csv") {
        metadata <- read.csv(metadata, header=TRUE)
    } else if (ext == "txt" || ext == "tsv") {
        metadata <- read.table(metadata, header=TRUE, sep="\t")
    } else {
        stop("Metadata file must be in .csv, .txt or .tsv format.")
    }

    return(metadata)
}

create_dds <- function(counts, metadata, sample_col, design, TE, level) {
    if (TE == TRUE && level == "subfamily") {
        counts$gene_id <- gsub("_dup[0-9]*", "", counts$gene_id)
        counts <- as.data.frame(counts %>% select(-ids) %>% group_by(gene_id) %>% summarize_all(sum))
    }
    dup <- (counts %>% group_by(gene_id) %>% summarize(n=n()) %>% filter(n>1))$gene_id
    counts <- counts[!(counts$gene_id %in% dup), ]
    counts$gene_id <- gsub("_(HS|MM)$", "", counts$gene_id)
    rownames(counts) <- counts$gene_id
    counts <- subset(counts, select=metadata[[sample_col]])
    counts <- round(as.matrix(counts))


    dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design = as.formula(paste("~", design)))

    dds <- DESeq(dds, parallel=FALSE)
    return(dds)
}

sizeFactors_plot <- function(dds, metadata, sample_col, condition_col, TE, level, output_dir) {
    size_factors <- setNames(as.data.frame(dds$sizeFactor), c("sizeFactor"))
    size_factors$sample <- rownames(size_factors)
    size_factors <- left_join(size_factors, metadata, by=c("sample" = sample_col))
    
    if (TE == TRUE) {
        TE <- "_TE"
        level <- paste0(level, "_")
    } else {
        TE <- ""
        level <- ""
    }

    filename <- file.path(output_dir, paste0(level, "size_factors", TE,".pdf"))
    plot <- ggplot(size_factors, aes(x=sizeFactor, y=sample, color=.data[[condition_col]])) +
                geom_bar(stat="identity", fill="white")+
                theme_minimal()+
                labs(color="Conditions", title="Size factors of all samples", y="Samples")+
                scale_color_brewer(palette="Set2")
    ggsave(filename, plot)
}

get_palette <- function(data, condition_col) {
    length <- length(unique(data[[condition_col]]))

    if (length <3){
        palette <- brewer.pal(3, "Set2")[1:length]
    } else if (length >=3 && length <= 8) {
        palette <- brewer.set2(length)
    } else if (length > 8 && length <= 12) {
        palette <- brewer.paired(length)
    } else if (length > 12 && length <= 26) {
        palette <- unname(alphabet2(length))
    } else if (length > 26 && length <= 32) {
        palette <- glasbey(length)
    } else {
        stop("More than 32 colors to plot, no available color palette.")
    }
    return(palette)
}

PCA_plot <- function(dds, metadata, topN=1000, TE, level, sample_col, condition_col, output_dir) {
    rld <- vst(dds)
	rv <- rowVars(assay(rld), useNames=TRUE)
	select <- head(order(rv, decreasing = TRUE),n=topN)
	pca <- prcomp(t(assay(rld)[select, ]))
    
    if (TE == TRUE) {
        TE <- "_TE"
        level <- paste0(level, "_")
    } else {
        TE <- ""
        level <- ""
    }

    filename <- file.path(output_dir, paste0(level, "pca_var", TE, ".pdf"))
	pdf(filename)
	pca.var <- pca$sdev^2
	pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
	barplot(pca.var.per, main="Variance explained by components",
				xlab="Principal Component", ylab="Percent Variation")
	dev.off()
	
    pca.data <- data.frame(Sample=rownames(pca$x),
	X=pca$x[, 1],
	Y=pca$x[, 2])
	pca.data <- left_join(pca.data, metadata, by=c("Sample" = sample_col))
    
    palette <- get_palette(pca.data, condition_col)

    filename <- file.path(output_dir, paste0(level, "pca", TE, ".pdf"))
	pdf(filename)
	print(ggplot(data=pca.data, aes(x=X, y=Y, label=Sample, color=.data[[condition_col]])) +
                geom_text(size=3, vjust=2.5) +
                geom_point(size=3) +
                labs(color="Conditions", size=10, label="Sample") +
                xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
                ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
                theme_bw() +
                theme(aspect.ratio=1) +
                scale_color_manual(values=palette)+
                ggtitle(paste0("PCA analysis (", topN, " most variable genes)")))
	dev.off()
}

clustering_plot <- function(dds, TE, level, output_dir) {
    rld <- vst(dds)
    dist2 <- as.dist(1-cor(assay(rld), method="pearson"))
    clust <- hclust(dist2, method="complete")

    if (TE == TRUE) {
        TE <- "_TE"
        level <- paste0(level, "_")
    } else {
        TE <- ""
        level <- ""
    }


    filename <- file.path(output_dir, paste0(level, "clustering_dendrograms", TE,".pdf"))
    pdf(filename)
    plot(clust, cex=0.8, xlab="Samples")
    dev.off()
    
    filename <- file.path(output_dir, paste0(level, "pheatmap_clustering_pearson", TE, ".pdf"))
    pdf(filename)
    pheatmap(as.matrix(dist2), clustering_distance_rows=dist2, clustering_distance_cols=dist2,
            color=hcl.colors(100, "BluYl"), main="Sample clustering (1 - Pearson correlation)")
    dev.off()
}

save_normalized_counts <- function(dds, output_dir) {
    counts <- as.data.frame(log2(counts(dds, normalized=TRUE)+1))
    counts$gene_id <- rownames(counts)
    filename <- file.path(output_dir, "log_normalized_counts.csv")
    write.csv(counts, file=filename, row.names=FALSE, quote=FALSE)
}

volcano_plot <- function(data, analysis_name, logFC, use_padj, padj, output_dir) {
    if (use_padj == TRUE) {
        pval <- "padj"
    } else {
        pval <- "pvalue"
    }
    
    data <- data[!is.na(data$log2FoldChange)&!is.na(data[[pval]]), ]
    perc <- as.numeric(padj) * 100

    filename <- file.path(output_dir, paste0(analysis_name, "_volcano_", pval, perc, "perc_logFC-", logFC, ".pdf"))
    plot <- ggplot(data, aes(x=log2FoldChange, y=-log10(.data[[pval]]), color=deg))+
                    geom_point(size=1)+
                    theme_bw()+
                    theme(panel.grid.minor=element_blank(),
                        panel.grid.major=element_line(color="#F2F2F2"),
                        axis.title=element_text(size=15),
                        axis.text=element_text(size=14),
                        legend.title=element_text(size=14),
                        legend.text=element_text(size=13))+
                    geom_hline(yintercept=c(-log10(padj)), linetype="dashed", color="black")+
                    geom_vline(xintercept=c(-logFC, logFC), linetype="dashed", color="black")+
                    labs(color="Differential expression")+
                    scale_color_manual(values=c("Overexpressed"="#F46969",
                                                "Underexpressed"="#6868FF",
                                                "Not significant"="#6C6C6C"))
    ggsave(filename, plot, height=4, width=7)
}

identify_DEGs <- function(data, logFC, use_padj, padj, output_dir) {
    if (use_padj == TRUE) {
        data$deg <- ifelse(data$log2FoldChange>=logFC&data$padj<=padj, "Overexpressed",
                    ifelse(data$log2FoldChange<=(-logFC)&data$padj<=padj, "Underexpressed",
                    "Not significant"))
    } else {
        data$deg <- ifelse(data$log2FoldChange>=logFC&data$pvalue<=padj, "Overexpressed",
                    ifelse(data$log2FoldChange<=(-logFC)&data$pvalue<=padj, "Underexpressed",
                    "Not significant"))
    }
    return(data)
}

contrast <- function(dds, experiment, control, analysis_name, condition_col, logFC, use_padj, padj, plot_volcano, output_dir) {
    data <- as.data.frame(results(dds, cooksCutoff=FALSE,
                    contrast=c(condition_col, experiment, control), alpha=0.05))
    cols <- colnames(data)
    data$gene_id <- rownames(data)
    data <- subset(data, select=c("gene_id", cols))

    filename <- file.path(output_dir, paste0(analysis_name, "_unfiltered.txt"))
    write.table(data, file=filename, row.names = FALSE, sep = "\t", quote = FALSE)

    data <- data[(!is.na(data$log2FoldChange)&!is.na(data$pvalue)&!is.na(data$gene_id)), ]
    filename <- file.path(output_dir, paste0(analysis_name, "_sorted_logFC_geneID.txt"))
    write.table(data[order(-(data$log2FoldChange), data$pvalue), ], file=filename, row.names=FALSE, sep="\t", quote=FALSE)

    data <- identify_DEGs(data, logFC, use_padj, padj, output_dir)

    if (plot_volcano == TRUE) {
        volcano_plot(data, analysis_name, logFC, use_padj, padj, output_dir)
    }

    data <- data %>% filter(deg != "Not significant")
    if (use_padj == TRUE) {
        pval <- "padj"
    } else {
        pval <- "pvalue"
    }

    data <- subset(data, select=c("gene_id", "deg", "log2FoldChange", "pvalue", "padj"))
    data <- data[order(-data$log2FoldChange, data[[pval]]),]
    perc <- as.numeric(padj) * 100
    filename <- file.path(output_dir, paste0(analysis_name, "_sorted_geneID_", pval, perc, "%_logFC-", logFC, ".txt"))
    write.table(data, file=filename, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

gsea <- function(args, analysis_name, output_dir) {
    script <- "/u/gironnea/my_libraries/GSEA_analysis.R"
    args <- c("--analysis_name", analysis_name,
              "--output_dir", output_dir,
              "--file", file.path(output_dir, paste0(analysis_name, "_unfiltered.txt")),
              "--plot", args$plot_GSEA,
              "--pathway", args$pathway)
    run("Rscript", c(script, args))
}

write_summary <- function(args) {
    lines <- c(paste0("counts: ", args$counts),
               paste0("metadata: ", args$metadata),
               paste0("sample_col: ", args$sample_col),
               paste0("condition_col: ", args$condition_col),
               paste0("design: ", args$design),
               paste0("deg analysis: ", args$deg_analysis),
               paste0("comparisons: ", args$comparisons),
               paste0("logFC: ", args$logFC),
               paste0("use_padj: ", args$use_padj),
               paste0("padj/pvalue: ", args$padj),
               paste0("TE: ", args$TE),
               paste0("level: ", args$level),
               paste0("run_GSEA: ", args$run_GSEA),
               paste0("pathway: ", args$pathway),
               paste0("output_dir: ", args$output_dir))
    date <- format(Sys.time(), "%Y%m%d_%H%M")
    writeLines(lines, file.path(args$output_dir, paste0("run_summary_", date, ".txt")))
}

adjust_args <- function(args) {
    if (args$run_GSEA == FALSE) {
        args$plot_GSEA <- FALSE
    }

    if (args$deg_analysis == FALSE) {
        args$plot_volcano <- FALSE
    }

    if (grepl("/", args$output_dir) == FALSE) {
        args$output_dir <- paste0("./", args$output_dir)
    }

    return(args)
}

main <- function(args) {
    create_output_dir(args$output_dir)
    counts <- read_counts(args$counts)
    metadata <- read_metadata(args$metadata)
    if (!is.na(args$comparisons)) {
        comparisons <- read.csv(args$comparisons, header=FALSE)
    }

    dds <- create_dds(counts, metadata, args$sample_col, args$design, args$TE, args$level)
    if (args$plot_PCA == TRUE) {
        PCA_plot(dds, metadata, 1000, args$TE, args$level, args$sample_col, args$condition_col, args$output_dir)
    }

    if (args$plot_sizeFactors == TRUE) {
        sizeFactors_plot(dds, metadata, args$sample_col, args$condition_col, args$TE, args$level, args$output_dir)
    }

    if (args$plot_clustering == TRUE) {
        clustering_plot(dds, args$TE, args$level, args$output_dir)
    }

    if (args$save_counts == TRUE) {
        save_normalized_counts(dds, args$output_dir)
    }

    if (args$deg_analysis == TRUE) {

        for (i in seq_len(dim(comparisons)[1])) {
            experiment <- comparisons$V1[i]
            control <- comparisons$V2[i]
            analysis_name <- paste(experiment, "vs", control, sep="_")
            output_dir <- file.path(args$output_dir, analysis_name)
            create_output_dir(output_dir)

            contrast(dds, experiment, control, analysis_name, args$condition_col, args$logFC, args$use_padj, args$padj,
                    args$plot_volcano, output_dir)
            
            if (args$TE == FALSE && args$run_GSEA == TRUE) {
                gsea(args, analysis_name, output_dir)
            }
        }
    }

    write_summary(args)
}



parser <- initialize_parser()
args <- parse_args(parser)
verifications(args)
args <- adjust_args(args)
main(args)