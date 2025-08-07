set -e
ulimit -n 65536
JOBNAME=local_job-$(date +%Y%m%d%H%M%S)-$RANDOM

ml star/2.7.11a
ml samtools/1.19

# Default values
seq=SE
samples=()
pre=results
output_dir="."
threads=32


help() {
    echo "Usage: star.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_R1 [path_to_R2]>    File(s) to analyze"
    echo "--pre <sample_prefix>                 Name to give to the output file         [default:results]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<./[star]>]"
    echo "--seq <PE/SE>                         Sequencing type                         [default:SE]"
    echo "--threads                             Number of threads                       [default:32]"
    echo "--index                               Path to index"
    echo "--help                                Show this help"
    echo ""
    echo ""
    echo "NB. If two files are given (eg. for PE samples), --sample needs to be specified twice, otherwise an error will be thrown"
    exit 0
}


throw_error() {
    echo "Error: ${1}" >&2 
    help
}

# Modify default values
if [[ $# == 0 ]]
then
    throw_error "No samples given. Exiting"
    exit 1
else
    options=$(getopt -o "" --long seq:,sample:,pre:,output_dir:,threads:,index:,help -n "$0" -- "$@")
    eval set -- "$options"


    # Verify getopt output code
    if [[ $? != 0 ]]
    then
        throw_error "Execution error with getopt"
        exit 1
    
    # Parse arguments
    else
        while true; do
            case "$1" in
                --sample)
                    samples+=("$2"); shift; shift
                    ;;
                --pre)
                    pre="$2"; shift; shift
                    ;;
                --seq)
                    seq="$2"; shift; shift
                    ;;
                --output_dir)
                    output_dir="$2"; shift; shift
                    ;;
                --threads)
                    threads="$2"; shift; shift
                    ;;
                --index)
                    index="$2"; shift; shift
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

options_summary() {
    echo "Script: star.sh"
    echo "Sample(s): ${samples[@]}"
    echo "Prefix: $pre"
    echo "Output directory: $output_dir"
    echo "Sequencing: $seq"
    echo "Threads: $threads"
    echo "Index: $index"
}

recap_message() {
    options_summary > $output_dir/${pre}_run_summary_$(date +%Y%m%d_%H%M).out
    touch $output_dir/$JOBNAME-${HOSTNAME%%.*}.done
}

verifications() {
    if [[ "${#samples[@]}" -eq 0 ]]
    then
        throw_error "No samples given. Exiting"
        exit 1
    fi
}

initialize_variables() {
    if [[ $output_dir == */ ]]; then
        output_dir=${output_dir%/}
    fi

    output_dir=$output_dir/star/$pre

    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi


    touch $output_dir/$JOBNAME-${HOSTNAME%%.*}.running
}


verifications
initialize_variables

STAR_options=(
        --runMode alignReads \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outStd Log \
        --runThreadN $threads \
        --outReadsUnmapped Fastx \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --outMultimapperOrder Random \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --limitBAMsortRAM 31532137230 \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN $threads \
        --outSAMunmapped Within \
        --quantMode GeneCounts \
        --outFileNamePrefix $output_dir/${pre}_ \
        --outSAMattributes NH HI AS NM MD
)

if [[ $seq == "PE" ]]; then
        STAR ${STAR_options[@]} \
            --genomeDir $genome \
            --readFilesIn ${samples[0]} ${samples[1]}
elif [[ $seq == "SE" ]]; then
        STAR ${STAR_options[@]} \
            --genomeDir $genome \
            --readFilesIn ${samples[0]}
fi

recap_message

echo $output_dir