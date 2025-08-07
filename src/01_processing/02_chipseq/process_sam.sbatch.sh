#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=150gb

BASE_DIR=./src/01_processing/CODES

set -e

ml samtools/1.19

#Default values
specie=hs
output_dir="."
seq="PE"
PE=true

help() {
    echo "Usage: process_sam.sbatch.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_sam>                File to filter"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<.>]"
    echo "--specie <hs/hsmm/hsmmdm/hsdm/hsec>   Species for reference genome selection  [default:hs]"
    echo "--seq <PE/SE>                         Sequencing type                         [default:PE]"
    echo "--help                                Show this help"
    exit 0
}

throw_error() {
    echo "Error: ${1}" >&2
    help
}

options_summary() {
    echo "Script: process_sam.sbatch.sh"
    echo "Sample: $sample"
    echo "Output directory: $output_dir"
    echo "Specie: $specie"
    echo "Sequencing: $seq"
}

recap_message() {
    options_summary > $output_dir/run_summary_$(date +%Y%m%d_%H%M).out
    touch $output_dir/$SLURM_JOB_ID-${HOSTNAME%%.*}.done
}

# Modify arguments
if [[ $# -eq 0 ]]
then
    throw_error "No sample given. Exiting"
    exit 1

else
    options=$(getopt -o "" --long specie:,sample:,output_dir:,seq:,help -n "$0" -- "$@")
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
                --specie)
                    specie="$2"; shift; shift
                    ;;
                --output_dir)
                    output_dir="$2"; shift; shift
                    ;;
                --seq)
                    seq="$2"; shift; shift
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

verifications() {
    if [[ $output_dir == */$* ]]; then
        output_dir=${output_dir%/*}
    fi

    case "$specie" in
        hs|hsmm|hsmmdm|hsdm|hsec) ;;
        *)
            throw_error "Specie not available"
            exit 1
    esac

    case "$seq" in
        SE|PE) ;;
        *)
            throw_error "Invalid sequencing type"
            exit 1
    esac

    if [[ ! $sample == *sam* ]]; then
        throw_error "Files must be in .sam format"
        exit 1
    fi
}

initialize_variables() {
    if [[ $output_dir == */ ]]; then
        output_dir=${output_dir%/}
    fi

    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi

    if [[ $seq == "SE" ]]; then
        PE=false
    fi

    touch $output_dir/$SLURM_JOB_ID-${HOSTNAME%%.*}.running
}

verifications
initialize_variables



echo "Filtering sam..." >&2 
julia $BASE_DIR/filter_sam.jl --file $sample --output_dir $output_dir --specie $specie --PE $PE

echo "Converting sam to bam..." >&2 
samtools view -u -b $output_dir/hs.sam -o $output_dir/hs.bam

echo "Sorting bam..." >&2 
samtools sort $output_dir/hs.bam > $output_dir/hs_sorted.bam

echo "Indexing bam..." >&2 
samtools index $output_dir/hs_sorted.bam -o $output_dir/hs.bam.bai


echo "File processed successfully"

recap_message