BASE_DIR=./src/01_processing/CODES

set -e

# Default values
samples=()
run_alignment=TRUE
run_trimming=TRUE
seq=SE
pre=results
output_dir="."
threads=32
TE=FALSE

help() {
    echo "Usage: RNAseq_pipeline.sbatch.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_R1 [path_to_R2]>    File(s) to analyze"
    echo "--pre                                 Name to give to the output file         [default:results]"
    echo "--alignment <TRUE/FALSE>              Align samples to transcriptome          [default:TRUE]"
    echo "--trimming <TRUE/FALSE>               Trim sequences                          [default:TRUE]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<./>]"
    echo "--seq <PE/SE>                         Sequencing type                         [default:SE]"
    echo "--threads                             Number of threads for alignment         [default:32]"
    echo "--index                               Path to index"
    echo "--TE <TRUE/FALSE>                     If TRUE, run TE alignment               [default:FALSE]"
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

# Verification of parameters
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

run_trimmomatic_function() {
    echo "Trimming sequences..."

    trimmed=($(bash $BASE_DIR/trimmomatic.sh ${samples[@]} --seq $seq --pre $pre --output_dir $output_dir))


    if [[ $? -ne 0 ]]; then
        throw_error "There was a problem with the processing"
        exit 1
    else
        echo "Sample(s) trimmed successfully. Trimmed sample(s) can be found here: "
    fi
    echo "${trimmed[@]}"
}


run_star_function() {
    echo "Aligning samples to reference genome..."

    done=$(bash $BASE_DIR/star.sh ${samples[@]} --pre $pre --seq $seq --output_dir $output_dir --threads $threads --index $index)
            
    if [[ $? -eq 0 ]]; then
        echo "Sample(s) aligned successfully. Aligned sample(s) can be found here: $output_dir/star/$pre"
    else
        throw_error "There was a problem with the alignment"
        exit 1
    fi
}

run_salmon_function() {
    echo "Aligning samples to reference genome..."

    done=$(bash $BASE_DIR/salmon.sh ${samples[@]} --pre $pre --seq $seq --output_dir $output_dir --threads $threads)
            
    if [[ $? -eq 0 ]]; then
        echo "Sample(s) aligned successfully. Aligned sample(s) can be found here: $output_dir/salmon/$pre"
    else
        throw_error "There was a problem with the alignment"
        exit 1
    fi

}

# Modify default values
if [[ $# -eq 0 ]]
then
    throw_error "No samples given. Exiting"
    exit 1

else
    options=$(getopt -o "" --long alignment:,trimming:,seq:,threads:,TE:,sample:,pre:,output_dir:,index:,help -n "$0" -- "$@")
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
                --alignment)
                    run_alignment="$2"; shift; shift
                    ;;
                --trimming)
                    run_trimming="$2"; shift; shift
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
                --TE)
                    TE="$2"; shift; shift
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

if [[ $output_dir == */ ]]; then
    output_dir=${output_dir%/}
fi

# Prepare options (expanded)
samples=("${samples[@]/#/--sample }")

if [[ $run_trimming == TRUE ]]
then
    run_trimmomatic_function

    if [[ $run_alignment == TRUE ]]; then
        samples=("${trimmed[@]/#/--sample }")
    
        if [[ $TE == FALSE ]]; then
            run_star_function
        else
            run_salmon_function
        fi
    fi
elif [[ $run_alignment == TRUE ]]; then
    if [[$TE == FALSE ]]; then
        run_star_function
    else
        run_salmon_function
    fi
fi