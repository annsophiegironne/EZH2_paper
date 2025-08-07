set -e

output_dir=salmon_index
specie=hs
kmer=31

help() {
    echo "Usage: index_salmon.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--fasta                               Fasta file"
    echo "--kmer                                Kmer size                               [default:31]"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<salmon_index>]"
    echo "--specie <hs/hsmm>                    Species for reference genome selection  [default:hs]"
    echo "--salmon                              Path to Salmon exe"
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
    options=$(getopt -o "" --long fasta:,output_dir:,specie:,salmon:,kmer:,help -n "$0" -- "$@")
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
                --kmer)
                    kmer="$2"; shift; shift
                    ;;
                --specie)
                    specie="$2"; shift; shift
                    ;;
                --salmon)
                    salmon="$2"; shift; shift
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
    echo "Script: index_salmon.sh"
    echo "fasta: $fasta"
    echo "Output directory: $output_dir"
    echo "Specie: $specie"
    echo "Salmon: $salmon"
}

recap_message() {
    options_summary > $output_dir/run_summary_$(date +%Y%m%d_%H%M).out
}

verifications() {
    case "$specie" in
        hs|hsmm) ;;
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
        hs)
            dir=Homo_sapiens
            ;;
        hsmm)
            dir=Homo_sapiens-Mus_musculus
            ;;
        *)
            throw_error "Unrecognized specie"
            exit 1
    esac
    version=$($salmon --version | cut -d ' ' -f 2)
    output_dir=$output_dir/$dir/salmon.v$version.$kmer-mer

    if [[ ! -d $output_dir ]]; then
        mkdir -p $output_dir
    fi

    if [ ! -z $(ls $output_dir) ]; then
        throw_error "Output directory $output_dir is not empty"
        exit 1
    fi
}


verifications
initialize_variables

$salmon index -t $fasta --i $output_dir -k $kmer --keepDuplicates --no-clip

recap_message
