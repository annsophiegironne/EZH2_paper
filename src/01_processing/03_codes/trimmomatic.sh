set -e

ml trimmomatic/0.39
JOBNAME=$SLURM_JOB_ID

# Default values
seq=PE
samples=()
pre=results
output_dir="."


help() {
    echo "Usage: trimmomatic.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--sample <path_to_R1 [path_to_R2]>    File(s) to analyze"
    echo "--pre <sample_prefix>                 Name to give to the output file         [default:results]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<./[trimmomatic]>]"
    echo "--seq <PE/SE>                         Sequencing type                         [default:SE]"
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

options_summary() {
    echo "Script: trimmomatic.sh"
    echo "Sample(s): ${samples[@]}"
    echo "Prefix: $pre"
    echo "Output directory: $output_dir"
    echo "Sequencing: $seq"
}

recap_message() {
    options_summary > $output_dir/${pre}_run_summary_$(date +%Y%m%d_%H%M).out
    touch $output_dir/$JOBNAME.done
}

error_message() {
    touch $output_dir/${pre}_run_summary_$(date +%Y%m%d_%H%M)_ERROR.out
}

# Modify default values of specified arguments
if [[ $# -eq 0 ]]
then
    throw_error "No samples given or output_dir specified. Exiting"
    exit 1
else
    options=$(getopt -o "" --long seq:,sample:,pre:,output_dir:,help -n "$0" -- "$@")
    eval set -- "$options"


    # Verify getopt output code
    if [[ $? != 0 ]]
    then
        throw_error "Execution error with getopt"
    
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

    output_dir=$output_dir/trimmomatic/$pre

    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi

    touch $output_dir/$JOBNAME.running
}

verifications
initialize_variables

touch $output_dir/$JOBNAME.running

if [[ $seq == "PE" ]]
then
    time java -Xms4g -Xmx4g -jar $TRIMMOMATIC_JAR PE \
                            -threads 32 -phred33 \
                            ${samples[0]} \
                            ${samples[1]} \
                            $output_dir/${pre}_R1_paired.fastq.gz \
                            $output_dir/${pre}_R1_unpaired.fastq.gz \
                            $output_dir/${pre}_R2_paired.fastq.gz \
                            $output_dir/${pre}_R2_unpaired.fastq.gz \
                            ILLUMINACLIP:$TRIMM_ADAPTER_DIR/TruSeq3-PE-2.fa:2:30:10:8:true \
                            LEADING:20 TRAILING:20 MINLEN:20 \
                            && echo $output_dir/${pre}_R1_paired.fastq.gz $output_dir/${pre}_R2_paired.fastq.gz
elif [[ $seq == "SE" ]]
then
    time java -Xms4g -Xmx4g -jar $TRIMMOMATIC_JAR SE \
                            -threads 32 -phred33 \
                            ${samples[0]} \
                            $output_dir/${pre}_trimmed.fastq.gz \
                            ILLUMINACLIP:$TRIMM_ADAPTER_DIR/TruSeq3-SE.fa:2:30:10 \
                            LEADING:20 TRAILING:20 MINLEN:20 \
                            && echo $output_dir/${pre}_trimmed.fastq.gz
fi

if [[ $? != 0 ]]; then
    error_message
    exit 1
else
    recap_message
fi