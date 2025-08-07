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
        "--specie"
            help="Species of reference genome"
            arg_type = String
            default = "hs"
        "--output_dir"
            help="Directory where to write output files. By default, output files will be written in current directory."
            arg_type = String
            default = "."
        "--PE"
            help="Whether reads are paired-end (PE) or not."
            arg_type = Bool
            default = true
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println(parsed_args)

function verifications(specie, file, output_dir)
    if !(specie in ["hs", "hsmm", "hsmmdm", "hsdm", "hsec"])
        throw("Specie not recognized. Available options are: hs, hsmm, hsmmdm, hsdm, hsec")
    end

    if !isfile(file)
        throw("Input file does not exit. Check path?")
    end

    if !endswith(file, ".sam")
        throw("Input file must be in sam format")
    end

    if !isfile(output_dir)
        mkpath(output_dir)
    end
end




function split_string(string, jump)
    strings=Vector{String}()

    for i in 1:jump:length(string)
        push!(strings, string[i:min(i+jump-1, length(string))])
    end

    return strings
end


function filter_alignment(line, MAPQ, sep, PE, output_paths, writing_files)
    if !occursin(r"[X|S]A:Z", line)
        line = split(line, sep, keepempty=false)
        if parse(Int16, line[5]) > MAPQ
            flag = parse(Int16, line[2])
            if flag < 256 && flag != 121 && flag != 137
                if (PE && line[7] == "=") || (PE == false)
                    line[3], specie = split_with_substring(line[3], "_")
                    index = findfirst(x -> occursin(lowercase(specie), x), output_paths)
                    println(writing_files[index], join(line, sep))
                end
            end
        end
    end
end

function split_with_substring(my_string, sep)
    index = findfirst(sep, my_string)[1]
    return SubString(my_string, 1, index-1), SubString(my_string, index+1, length(my_string))
end


function process_alignments(input_file, sep, species, PE, output_dir)
    output_paths = [joinpath(output_dir, "$(x).sam") for x in species]
    writing_files = [open(f, "w") for f in output_paths]

    open(String(input_file), "r") do file
        for s in eachline(file)
            if startswith(s, "@SQ")
                s = split(s, sep, keepempty=false)
                s[2], specie = split_with_substring(s[2], "_")
                index = findfirst(x -> occursin(lowercase(specie), x), output_paths)
                if index != nothing
                    println(writing_files[index], join(s, sep))
                else
                    throw("Unexpected specie chromosome found in file. Check settings")
                end
            elseif startswith(s, "@PG")
                for writing_file in writing_files
                    println(writing_file, s)
                end
            else
                # First non-header line
                filter_alignment(s, 30, '\t', PE, output_paths, writing_files)
                break
            end
        end

        for s in eachline(file)
            filter_alignment(s, 30, '\t', PE, output_paths, writing_files)
        end
    end

    for output_file in writing_files
        close(output_file)
    end
end

function check_output_dir(output_dir)
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    return output_dir
end

function main(file, specie, PE, output_dir)
    output_dir = verifications(specie, file, output_dir)
    output_dir = check_output_dir(output_dir)
    species = split_string(specie, 2)
    @time process_alignments(file, '\t', species, PE, output_dir)
end

main(parsed_args["file"], parsed_args["specie"], parsed_args["PE"], parsed_args["output_dir"])