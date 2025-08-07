using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--files"
            help = "Mean will be computed on input bw files"
            arg_type = String
            nargs = '*'
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

function verifications(files)
    for file in files
        if !endswith(file, "bw") && !endswith(file, "bigwig")
            throw("Input file must be in bw format (.bw, .bigwig)")
        end
    end
end

function create_output_dir(file)
    output_dir = dirname(file)
    if !isdir(output_dir)
        mkpath(output_dir)
    end
end

function main(files, convert_to_bw, output_file)
    verifications(files)
    create_output_dir(output_file)
    if isfile(output_file)
        run(`rm $output_file`)
    end
    run( `wiggletools/bin/wiggletools write_bg $output_file mean $files`)
    if convert_to_bw == true
        run(pipeline(` sort -k1,1 -k2,2n $output_file`, stdout= output_file*".tmp"))
        bw = replace(output_file, r"\.bdg" => ".bw")
        run(` ucsc-v396/bin/bedGraphToBigWig $output_file.tmp src/01_processing/03_codes/annotations/Homo_sapiens.GRCh38.chrom_sizes $bw `)
        run(` rm $output_file.tmp `)
    end 
end


@time main(parsed_args["files"], parsed_args["convert_to_bw"], parsed_args["output_file"])