set -e

# Default parameters
output_dir=star_index
bp=99
specie=hs
gencode=40
version=$(STAR --version)

help() {
    echo "Usage: index_star.sh [OPTIONS]"
    echo ""
    echo "Options:"
    echo "--gtf                                 GTF annotation file"
    echo "--fasta                               Fasta file"
    echo "--output_dir <path_to_dir>            Directory where to write output files   [default:<star_index>]"
    echo "--bp <int>                            Overhang (length of reads - 1)          [default:99]"
    echo "--specie <hs/hsmm>                    Species for reference genome selection  [default:hs]"
    echo "--gencode                             Version of Gencode used                 [default:40]"
    echo "--help                                Show this help"
    exit 0
}

throw_error() {
    echo "Error: ${1}" >&2 
    help
}


if [[ $# == 0 ]]
then
    throw_error "No GTF or fasta files given. Exiting"
    exit 1
else
    options=$(getopt -o "" --long gtf:,fasta:,output_dir:,bp:,specie:,gencode:,help -n "$0" -- "$@")
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
                --gtf)
                    gtf="$2"; shift; shift
                    ;;
                --fasta)
                    fasta="$2"; shift; shift
                    ;;
                --output_dir)
                    output_dir="$2"; shift; shift
                    ;;
                --bp)
                    bp="$2"; shift; shift
                    ;;
                --specie)
                    specie="$2"; shift; shift
                    ;;
                --gencode)
                    gencode="$2"; shift; shift
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
    echo "Script: index_star.sh"
    echo "GTF: $gtf"
    echo "fasta: $fasta"
    echo "Output directory: $output_dir"
    echo "Overhang: $bp"
    echo "Specie: $specie"
    echo "Gencode: $gencode"
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

    if [[ $fasta == "" ]]; then
        throw_error "No fasta file given"
    fi

    if [[ $fasta != *".fasta" && $fasta != *".fa" ]]; then
        throw_error "Unexpected fasta format"
    fi

    if [[ $gtf != *".gtf" ]]; then
        throw_error "Unexpected gtf format"
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

    output_dir=$output_dir/$dir/Gencode$gencode.STAR.v$version.sjdbOverhang$bp
    
    mkdir -p $output_dir

    if [[ ! -z $(ls $output_dir) ]]; then
        throw_error "Output directory $output_dir is not empty"
        exit 1
    fi
}


verifications
initialize_variables


STAR_options=(
    --runMode genomeGenerate --runThreadN 16 \
    --genomeFastaFiles $fasta \
    --genomeSAindexNbases 14 \
    --sjdbGTFfile $gtf \
    --limitGenomeGenerateRAM 50000000000000 \
    --sjdbOverhang $bp
)


STAR ${STAR_options[@]} --genomeDir $output_dir --outFileNamePrefix $output_dir

recap_message