set -e

ml bwa/0.7.17
JOBNAME=local_job-$(date +%Y%m%d%H%M%S)-$RANDOM


# Valeurs par défaut
seq=SE
samples=()
pre=results
output_dir="."
threads=32


help() {
    echo "Usage: bwa.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_R1 [path_to_R2]>    File(s) to analyze"
    echo "--pre <sample_prefix>                 Name to give to the output file         [default:results]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<./>]"
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
if [[ $# -eq 0 ]]
then
    throw_error "No samples given. Exiting"
    exit 1
else
    options=$(getopt -o "" --long seq:,sample:,pre:,output_dir:,threads:,index:,help -n "$0" -- "$@")
    eval set -- "$options"


    # Vérifie le code de sortie de getopt
    if [[ $? != 0 ]]
    then
        throw_error "Execution error with getopt"
        exit 1
    
    # Parse les arguments
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
    echo "Script: bwa.sh"
    echo "Sample(s): ${samples[@]}"
    echo "Prefix: $pre"
    echo "Output directory: $output_dir"
    echo "Sequencing: $seq"
    echo "Threads: $threads"
    echo "Index: $index"
}

recap_message() {
    options_summary > $output_dir/${pre}_run_summary_$(date +%Y%m%d_%H%M).out
    touch $output_dir/$JOBNAME.done
}


verifications() {
    if [[ "${#samples[@]}" -eq 0 ]]
    then
        throw_error "No samples given. Exiting"
        exit 1
    fi

    if [[ $seq == "PE" && ${#samples[@]} -ne 2 ]]; then
        throw_error "Two fastq files needed with PE samples"
        exit 1
    elif [[ $seq == "SE" && ${#samples[@]} -ne 1 ]]; then
        throw_error "One fastq file needed with SE samples"
        exit 1
    fi

    if [[ ! ${samples[0]} == *.fastq* ]]; then
        throw_error "Files must be in .fastq format (gzipped or not)"
        exit 1
    fi

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

    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi

    touch $output_dir/$JOBNAME.running
}

verifications
initialize_variables

if [[ $seq == "PE" ]]; then
    time bwa mem $genome -t $threads ${samples[0]} ${samples[1]} > $output_dir/${pre}.sam
elif [[ $seq == "SE" ]]; then
    time bwa mem $genome -t $threads ${samples[0]} > $output_dir/${pre}.sam
fi

recap_message

echo $output_dir/${pre}.sam