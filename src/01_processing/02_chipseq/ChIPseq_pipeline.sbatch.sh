BASE_DIR=.src/01_processing/CODES

set -e

# Default values
samples=()
run_BWA=TRUE
run_Trimmomatic=TRUE
seq=PE
pre=results
output_dir="."
threads=32


help() {
    echo "Usage: ChIPseq_pipeline.sbatch.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_R1 [path_to_R2]>    File(s) to analyze"
    echo "--pre <sample_prefix>                 Name to give to the output file         [default:results]"
    echo "--bwa <TRUE/FALSE>                    Align samples to genome                 [default:TRUE]"
    echo "--trimmomatic <TRUE/FALSE>            Trim sequences                          [default:TRUE]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<./>]"
    echo "--seq <PE/SE>                         Sequencing type                         [default:PE]"
    echo "--threads                             Number of threads for alignment         [default:32]"
    echo "--index                               Path to index"
    echo "--help                                Show this help"
    echo ""
    echo ""
    echo "NB. If two files are given (eg. for PE samples), --sample needs to be specified twice, otherwise an error will be thrown."
    echo "Fastq files are expected to be given in order for PE samples (R1 first, then R2), however this is not checked."
    exit 0
}

throw_error() {
    echo "Error: ${1}" >&2
    help
}

# Verify input parameters
verifications() {
    case "$seq" in
        SE|PE) ;;
        *)
            throw_error "Sequencing parameter not recognized"
            exit 1
    esac

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

    # Prepare options (expanded)
    samples=("${samples[@]/#/--sample }")
}

run_trimmomatic_function() {
    echo "Trimming sequences..." >&2

    trimmed=($( bash $BASE_DIR/trimmomatic.sh ${samples[@]} --seq $seq --pre $pre --output_dir $output_dir ))


    if [[ $? -ne 0 ]]; then
        throw_error "There was a problem with the processing"
        exit 1
    else
        echo "Sample(s) trimmed successfully. Trimmed sample(s) can be found here: $output_dir/trimmomatic/$pre" >&2
    fi
    # Return trimmed fastq paths
    echo ${trimmed[@]}
}


run_bwa_function() {
    echo "Aligning samples to reference genome..." >&2

    done=$(bash $BASE_DIR/bwa.sh ${samples[@]} --genome $genome --pre $pre --seq $seq --output_dir $output_dir/bwa --threads $threads)
            
    if [[ $? -eq 0 ]]; then
        echo "Sample(s) aligned successfully. Aligned sample(s) can be found here: $output_dir/bwa/$pre" >&2
    else
        throw_error "There was a problem with the alignment"
        exit 1
    fi
}

# Modify default alues
if [[ $# -eq 0 ]]
then
    throw_error "No samples given. Exiting"
    exit 1

else
    options=$(getopt -o "" --long bwa:,trimmomatic:,seq:,sample:,pre:,output_dir:,threads:,index:,help -n "$0" -- "$@")
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
                    samples+=("$2"); shift; shift
                    ;;
                --pre)
                    pre="$2"; shift; shift
                    ;;
                --bwa)
                    run_BWA="$2"; shift; shift
                    ;;
                --trimmomatic)
                    run_Trimmomatic="$2"; shift; shift
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
                    index="$2"; shift;shift
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



if [[ $run_Trimmomatic == TRUE ]]
then
    trimmed=($(run_trimmomatic_function))


    if [[ $run_BWA == TRUE ]]; then
        samples=("${trimmed[@]/#/--sample }")
        run_bwa_function
    fi
elif [[ $run_BWA == TRUE ]]; then
    run_bwa_function
fi