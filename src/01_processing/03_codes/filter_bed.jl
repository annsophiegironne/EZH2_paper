using ArgParse
using CSV
using DataFrames

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file"
            help="Alignment file in sam format to filter"
            arg_type = String
            required = true
        "--analysis_name"
            help="Name to give to the output file"
            arg_type = String
        "--output_dir"
            help="Directory where to write output files. By default, output files will be written in the same directory as the input file."
            arg_type = String
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println(parsed_args)


function verifications(file, analysis_name, output_dir)
    if output_dir == nothing || output_dir == ""
        if findlast('/', file) == nothing
            output_dir="./"
        else
            output_dir=SubString(file, 1, findlast('/', file))
        end
    end

    if !endswith(file, "bed") && !endswith(file, "Peak")
        throw("Input file must be in bed format (.bed, .narrowPeak or .broadPeak)")
    end

    if analysis_name == nothing || analysis_name == ""
        ind = findall('/', file)
        analysis_name = SubString(file, ind[end-1]+1, ind[end]-1)
    end

    return output_dir, analysis_name
end

function check_output_dir(output_dir)
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    return output_dir
end

function filter_bed(input, analysis_name, output_dir)
    chrs = Vector{String}(undef, 25)
    j = 1
    for i in [1:22; 'X'; 'Y'; 'M']
        chrs[j] = "chr$(i)"
        j+=1
    end

    output_file=joinpath(output_dir, "$(analysis_name).bed")
    writing_file=open(output_file, "w")

    open(input, "r") do file
        for s in eachline(file)
            ind = findall('\t', s)
            if SubString(s, 1, ind[1]-1) in chrs
                println(writing_file, SubString(s, 1, ind[3]-1))
            end
        end
    end

    close(writing_file)
end


function main(file, analysis_name, output_dir)
    output_dir, analysis_name = verifications(file, analysis_name, output_dir)
    output_dir = check_output_dir(output_dir)
    filter_bed(file, analysis_name, output_dir)
end


main(parsed_args["file"], parsed_args["analysis_name"], parsed_args["output_dir"])
