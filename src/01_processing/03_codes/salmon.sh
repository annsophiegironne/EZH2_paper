set -e

JOBNAME=local_job-$(date +%Y%m%d%H%M%S)-$RANDOM

# Default values
seq=SE
samples=()
pre=results
output_dir="."
threads=32

help() {
    echo "Usage: salmon.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_R1 [path_to_R2]>    File(s) to analyze"
    echo "--pre <sample_prefix>                 Name to give to the output file                [default:results]"
    echo "--output_dir <path_to_dir>            Directory where to write output files          [default:<./[salmon]>]"
    echo "--seq <PE/SE>                         Sequencing type                                [default:SE]"
    echo "--threads                             Number of threads                              [default:32]"
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

# Modifie les valeurs par défaut des arguments spécifiés
if [[ $# == 0 ]]
then
    throw_error "No samples given or output_dir specified. Exiting"
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
    echo "Script: salmon.sh"
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

    output_dir=$output_dir/salmon/$pre

    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi

    touch $output_dir/$JOBNAME-${HOSTNAME%%.*}.running
}


verifications
initialize_variables


if [[ $seq == "PE" ]]; then
        $salmon quant -i $index -l A -1 ${samples[0]} -2 ${samples[1]} -p $threads --softclip --recoverOrphans -o $output_dir
elif [[ $seq == "SE" ]]; then
        $salmon quant -i $index -l A -r ${samples[0]} -p $threads --softclip --recoverOrphans -o $output_dir
fi

recap_message

echo $output_dir/$pre