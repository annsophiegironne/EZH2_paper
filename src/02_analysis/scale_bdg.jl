using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file"
            help = "Bdg file to scale"
            arg_type = String
            required = true
        "--factor"
            help = "Scaling factor to apply"
            arg_type = Float16
            required = true
        "--convert_to_bw"
            help = "If true, scaled bdg file will also be outputted as a bw file"
            arg_type = Bool
            required = false
            default = true
        "--output"
            help = "Name of the output file. If directory not specified in name, output file will be written in current directory"
            default = "./results.bdg"
            dest_name = "output_file"
            arg_type = String
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println(parsed_args)


function verifications(file, factor)
    if !endswith(file, "bdg") && !endswith(file, "bedgraph")
        throw("Input file must be in bdg format (.bdg, .bedgraph)")
    end

    if factor <= 0
        throw("Normalization factor must be positive")
    end
end


function normalize_bdg(file, factor, output_file)
    open(file, "r") do input
        open(output_file, "w") do output
            for s in eachline(input)
                s = split(s, "\t")
                s[4] = string(parse(Float64, s[4]) / factor)
                s = join(s, "\t")
                println(output, s)
            end
        end
    end
end

function create_output_dir(file)
    output_dir = dirname(file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
end

function main(file, factor, output_file, convert_to_bw)
    verifications(file, factor)
    create_output_dir(output_file)
    normalize_bdg(file, factor, output_file)
    if convert_to_bw == true
        bw = replace(output_file, r"\.bdg" => ".bw")
        run(` ucsc-v396/bin/bedGraphToBigWig $output_file src/01_processing/03_codes/annotations/Homo_sapiens.GRCh38.chrom_sizes $bw `)
    end
end


main(parsed_args["file"], parsed_args["factor"], parsed_args["output_file"], parsed_args["convert_to_bw"])