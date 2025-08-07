library(ChIPseeker)
library(argparser)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

initialize_parser <- function() {
    parser <- arg_parser("Annotate peaks in BED file")
    parser <- add_argument(parser, "--file", type="character",
                        help="Path to BED file. This file should not have a header.")
    return(parser)
}


main <- function(args) {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene   
    peaks <- read.table(args$file, sep="\t", header=FALSE)
    peaks_GR <- GRanges(seqnames=peaks$V1, ranges=IRanges(peaks$V2, peaks$V3))
    
    anno <- annotatePeak(peaks_GR, TxDb=txdb, verbose=FALSE, annoDb="org.Hs.eg.db", columns=c("SYMBOL"))
    
    filename <- gsub(".bed", "_annotations.txt", args$file)
    write.table(anno, file=filename, sep="\t", quote=FALSE, row.names=FALSE)
}


parser <- initialize_parser()
args <- parse_args(parser)
main(args)