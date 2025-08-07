library(DiffBind)
library(dplyr)
library(tidyr)
library(DESeq2)
library(argparser)
library(tools)
library(stringr)

initialize_parser <- function() {
    parser <- arg_parser("Calculate differential binding analysis using DiffBind")
    parser <- add_argument(parser, "--metadata", type="character",
                          help="Path to metadata. This file should include all the samples that will be used for normalization")
    parser <- add_argument(parser, "--analysis_name", type="character", default="results",
                          help="Comparison identification, eg 'treatment_vs_control")
    parser <- add_argument(parser, "--groups", type="character",
                          help="A vector of two groups to be used for the binding analysis. 
                          Groups should be separated by a comma (e.g., 'treatment,control'). The last group will be used as control.")
    parser <- add_argument(parser, "--output_dir", type="character", default=".",
                          help="Path to the directory where results will be saved. If the directory does not exist, it will be created")
    parser <- add_argument(parser, "--minOverlap", type="numeric", default=1,
                          help="Minimum overlap between peaksets for merging and counting")
    parser <- add_argument(parser, "--importCounts", type="logical", default=FALSE,
                          help="If TRUE, import counts from dba object located at the path specified by --countsDir parameter")
    parser <- add_argument(parser, "--countsDir", type="character", default=NULL,
                          help="If --importCounts is set to TRUE, the dba count object will be found at this path")
    parser <- add_argument(parser, "--norm", type="character", default="DBA_NORM_LIB",
                          help="Normalization method to use.")
    parser <- add_argument(parser, "--factors", type="character", default="",
                          help="Path to normalization factors. This file should have no header, only have one column
                                and the normalization factors should be in the same order as the samples in the metadata.")
    return(parser)
}

verifications <- function(args) {
    if (is.na(args$metadata)) {
        stop("No metadata file given.")
    }

    if (is.na(args$groups)) {
        stop("No groups given.")
    }
}

MAplots <- function(DBA, analysis_name, run=FALSE, normalization=TRUE, output_dir) {
    if (run==TRUE) {
        if (normalization==TRUE) {
            filename <- file.path(output_dir, paste0(analysis_name, "_afterNormalization_MA.pdf"))
            pdf(filename)
            print(dba.plotMA(DBA, bNormalized=TRUE))
            dev.off()
        } else {
            filename <- file.path(output_dir, paste0(analysis_name, "_withoutNormalization_MA.pdf"))
            pdf(filename)
            print(dba.plotMA(DBA, bNormalized=FALSE))
            dev.off()
        }
    }
}


exportReport <- function(DBA, treatment, control, analysis_name, output_dir) {
    FCcutoff <- 1
    Pcutoff <- 0.05
    report <- dba.report(DBA, fold=0, th=1, bNormalized=TRUE, DataType=DBA_DATA_FRAME, bCalled=TRUE)
    report$Chr <- gsub("_HS", "", report$Chr)
    
    # Whole report
    filename <- file.path(output_dir, paste0(analysis_name, "_wholeReport.txt"))
    write.table(format(report, scientific=FALSE), file=filename, sep="\t", row.names=FALSE, quote=FALSE)
    
    #Specific to treatment
    filename <- file.path(output_dir, paste0(analysis_name, "_specific_to_", treatment, ".bed"))
    treatment <- report[(report$Fold > (FCcutoff)&report$"p-value"<Pcutoff), ]
    write.table(treatment[order(treatment$Fold, decreasing=TRUE), ], file=filename, sep="\t", row.names=FALSE, quote=FALSE)
    
    #Specific to control
    filename <- file.path(output_dir, paste0(analysis_name, "_specific_to_", control, ".bed"))
    control <- report[(report$Fold < (-FCcutoff)&report$"p-value"<Pcutoff), ]
    write.table(control[order(control$Fold, decreasing=TRUE), ], file=filename, sep="\t", row.names=FALSE, quote=FALSE)
    
    #Common
    common <- report[(report$Fold<FCcutoff&report$Fold>(-FCcutoff))|(report$"p-value">Pcutoff), ]
    filename <- file.path(output_dir, paste0(analysis_name, "_common.bed"))
    write.table(common[order(common$Fold, decreasing=TRUE), ], file=filename, sep="\t", row.names=FALSE, quote=FALSE)
    print("Reports exported")
}



contrast <- function(df, analysis_name, treatment, control, factors, normalization, mergeOverlap, countOverlap, importCounts, countsDir, output_dir) {
    if (importCounts==TRUE) {
        file <- paste0(analysis_name, "_counts")
        data <- dba.load(file=file, dir=countsDir)
    } else {
        data <- dba(sampleSheet=df, bRemoveRandom=TRUE, minOverlap=mergeOverlap, config=data.frame(fragmentSize=0, doBlacklist=FALSE, doGreylist=FALSE))
        data <- dba.count(data, bParallel=FALSE, summits=FALSE, minOverlap=countOverlap, bSubControl=TRUE, bScaleControl=TRUE)
        dba.save(data, file=paste0(analysis_name, "_counts"))
        print("Done counting and saved")
    }

    if (factors == "") {
        data <- dba.normalize(data, method=DBA_DESEQ2, normalize=get(normalization))

    } else {
        norm_factors <- read.table(factors, header=FALSE)$V1
        data <- dba.normalize(data, method=DBA_DESEQ2, normalize=norm_factors)

    }

    norm <- dba.normalize(data, method=DBA_DESEQ2, bRetrieve=TRUE)
    
    #Merge samples info with normalization factors
    norm_df <- dba.show(data)
    norm <- cbind(norm_df, norm)
    filename <- file.path(output_dir, paste0(analysis_name, "_normFactors.txt"))
    write.table(norm, file=filename, sep="\t", quote=FALSE, row.names=FALSE)
    print("Done normalizing")
    
    data <- dba.contrast(data, group1=data$masks[[treatment]], name1=treatment, group2=data$masks[[control]], name2=control)
    data <- dba.analyze(data, bParallel=TRUE, bBlacklist=FALSE, bGreylist=FALSE, method=DBA_DESEQ2)
    print("Done analyzing")
    
    exportReport(data, treatment, control, analysis_name, output_dir)
    
    MAplots(data, analysis_name, TRUE, TRUE, output_dir)
    MAplots(data, analysis_name, TRUE, FALSE, output_dir)
}

read_metadata <- function(file) {
    ext <- file_ext(file)

    if (ext == "txt") {
        data <- read.table(file, header=TRUE)
    } else if (ext == "csv") {
        data <- read.csv(file, header=TRUE)
    } else {
        stop("Metadata file should be txt or csv format.")
    }

    return(data)
}


create_output_dir <- function(output_dir) {
    if (!file.exists(output_dir)) {
        dir.create(output_dir, recursive=TRUE)
    }
}

write_summary <- function(args) {
    lines <- c(paste0("metadata: ", args$metadata),
               paste0("output_dir: ", args$output_dir),
               paste0("groups: ", args$groups),
               paste0("overlap: ", args$minOverlap),
               paste0("analysis_name: ", args$analysis_name),
               paste0("importCounts: ", args$importCounts),
               paste0("countsDir: ", args$countsDir),
               paste0("normalization: ", args$norm),
               paste0("factors: ", args$factors))
    date <- format(Sys.time(), "%Y%m%d_%H%M")
    writeLines(lines, file.path(args$output_dir, paste0("run_summary_", date, ".txt")))
}

main <- function(args) {
    df <- read_metadata(args$metadata)
    create_output_dir(args$output_dir)

    vector <- str_split_1(args$groups, ",")
    treatment <- vector[1]
    control <- vector[2]


    start <- Sys.time()
    contrast(df, args$analysis_name, treatment, control, args$factors, args$norm, args$minOverlap, args$minOverlap, args$importCounts, args$countsDir, args$output_dir)
    end <- Sys.time()
    time <- end-start
    print(time)

    write_summary(args)
}


parser <- initialize_parser()
args <- parse_args(parser)
verifications(args)
main(args)