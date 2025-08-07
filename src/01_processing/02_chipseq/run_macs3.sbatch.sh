#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=40gb

BASE_DIR=./src/01_processing/CODES

set -e

source env/bin/activate
chromsizes=$BASE_DIR/annotations/Homo_sapiens.GRCh38.chrom_sizes

ml ucsc/v396 

# Default values
pre="results"
output_dir="."
seq="PE"
control=""
fdr=0.01
broad_fdr=0.01

help() {
    echo "Usage: run_macs3.sbatch.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_bam>                File to analyze"
    echo "--control <path_to_bam>               Control file"
    echo "--mark <mark_or_target>               Histone mark or target of sample"
    echo "--pre <sample_prefix>                 Name to give to the output subdirectory [default:results]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<.>/[pre]]"
    echo "--seq [PE/SE]                         Sequencing type                         [default:PE]"
    echo "--fdr                                 FDR cut-off on peak calling             [default:0.01]"
    echo "--broad_fdr                           FDR cut-off on broad peak calling       [default:0.01]"
    echo "--help                                Show this help"
    echo ""
    echo ""
    echo ""
    echo "NB. FDR cut-off is applied to both narrow/broad peaks. However, broad peaks have a second peak calling cut-off,"
    echo "which can be set using --broad_fdr. See https://macs3-project.github.io/MACS/docs/callpeak.html for more info."
    exit 0
}

throw_error() {
    echo "Error: ${1}" >&2
    help
}

# Verify arguments
verifications() {
    if [[ ! $sample == *.bam$* ]]; then
        throw_error "Files must be in .bam format"
        exit 1
    fi

    if [[ $control != "" && $control != *.bam$* ]]; then
        throw_error "Control file must be in .bam format"
        exit 1
    fi
}

initialize_variables() {
    if [[ $output_dir == */ ]]; then
        output_dir=${output_dir%/}
    fi

    mkdir -p $output_dir/$pre 

    case "$mark" in
        H3K27me3|H3K9me3|H3K9me2|H3K36me3|H3K27me2|H3K36me2)
            width="broad"
            broad="--broad-cutoff $broad_fdr --broad"
            ;;
        *)
            width="narrow"
            broad=""
            ;;
    esac

    case "$seq" in
        SE)
            format="BAM" ;;
        PE)
            format="BAMPE" ;;
    esac

    if [[ $control != "" ]]; then
        control="${control/#/-c }"
    fi

    touch $output_dir/$pre/$SLURM_JOB_ID-${HOSTNAME%%.*}.running
}

options_summary() {
    control_sed=$(echo $control | sed 's/-c //g')
    echo "Script: run_macs3.sbatch.sh"
    echo "Sample: $sample"
    echo "Control: $control_sed"
    echo "Prefix: $pre"
    echo "Mark: $mark"
    echo "Output directory: $output_dir"
    echo "Sequencing: $seq"
    echo "FDR: $fdr"
    echo "Broad FDR: $broad_fdr"
}

recap_message() {
    options_summary > $output_dir/$pre/${pre}_run_summary_$(date +%Y%m%d_%H%M).out
    touch $output_dir/$pre/$SLURM_JOB_ID-${HOSTNAME%%.*}.done
}


# Modify arguments
if [[ $# -eq 0 || $# -eq 1 ]]
then
    throw_error "No samples or target given. Exiting"
    exit 1

else
    options=$(getopt -o "" --long sample:,pre:,mark:,control:,seq:,output_dir:,fdr:,broad_fdr:,help -n "$0" -- "$@")
    eval set -- "$options"


    # Verify getopt output code
    if [[ $? -ne 0 ]]
    then
        throw_error "Execution error with getopt"
        exit 1
    
    # Parse arguments
    else
        while true; do
            case "$1" in
                --sample)
                    sample="$2"; shift; shift
                    ;;
                --pre)
                    pre="$2"; shift; shift
                    ;;
                --mark)
                    mark="$2"; shift; shift
                    ;;
                --output_dir)
                    output_dir="$2"; shift; shift
                    ;;
                --control)
                    control="$2"; shift; shift
                    ;;
                --seq)
                    seq="$2"; shift; shift
                    ;;
                --fdr)
                    fdr=$2; shift; shift
                    ;;
                --broad_fdr)
                    broad_fdr=$2; shift; shift
                    ;;
                --help)
                    help ;;
                --)
                    shift 
                    ;;
                *)
                    break
            esac
        done
    fi
fi

verifications
initialize_variables

macs3 callpeak -t $sample $control --format $format -g hs -B $broad -q $fdr --outdir $output_dir/$pre -n $pre
julia $BASE_DIR/filter_bed.jl --file $output_dir/$pre/$pre*${width}Peak --analysis_name $pre --output_dir $output_dir/$pre

sort -k1,1 -k2,2n $output_dir/$pre/*treat_pileup.bdg > $output_dir/$pre/${pre}.sorted.bdg
bedGraphToBigWig $output_dir/$pre/${pre}.sorted.bdg $chromsizes $output_dir/$pre/$pre.bw

rm $output_dir/$pre/${pre}.sorted.bdg

recap_message
