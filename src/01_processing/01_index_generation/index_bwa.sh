set -e

# Default parameters
output_dir=genome_index
specie=hsdm
version=$(bwa 2>&1 | grep 'Version' | sed 's/.*: \(.*\)-.*/\1/g')

help() {
    echo "Usage: index_bwa.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--fasta                               Fasta file"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<genome_index>]"
    echo "--specie <hsdm/hsmmdm>                Species for reference genome selection  [default:hsdm]"
    echo "--help                                Show this help"
    exit 0
}

throw_error() {
    echo "Error: ${1}" >&2 
    help
}


if [[ $# == 0 ]]
then
    throw_error "No fasta file given. Exiting"
    exit 1
else
    options=$(getopt -o "" --long fasta:,output_dir:,specie:,help -n "$0" -- "$@")
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
                --fasta)
                    fasta="$2"; shift; shift
                    ;;
                --output_dir)
                    output_dir="$2"; shift; shift
                    ;;
                --specie)
                    specie="$2"; shift; shift
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
    echo "Script: index_bwa.sh"
    echo "fasta: $fasta"
    echo "Output directory: $output_dir"
    echo "Specie: $specie"
}

recap_message() {
    options_summary > $output_dir/run_summary_$(date +%Y%m%d_%H%M).out
}

verifications() {
    case "$specie" in
        hsdm|hsmmdm) ;;
        *)
            throw_error "Unexpected specie mentionned"
            exit 1
    esac

    if [[ $fasta != *".fasta" && $fasta != *".fa" ]]; then
        throw_error "Unexpected fasta format"
    fi

}

initialize_variables() {
    case "$specie" in
        hsdm)
            dir=Homo_sapiens-Drosophila_melanogaster
            ;;
        hsmmdm)
            dir=Homo_sapiens-Mus_musculus-Drosophila_melanogaster
            ;;
        *)
            throw_error "Unrecognized specie"
            exit 1
    esac

    output_dir=$output_dir/$dir/bwa.v$version
    pre=$output_dir/bwa.v$version

    mkdir -p $output_dir

    if [ ! -z $(ls $output_dir) ]; then
        throw_error "Output directory $output_dir is not empty"
        exit 1
    fi
}

copy_fasta() {
    new_fasta=$output_dir/bwa.v$version.fa
    cp $fasta $new_fasta
}

verifications
initialize_variables
copy_fasta


bwa index -p $pre $new_fasta

recap_message
