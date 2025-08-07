using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--file"
            help = "Path to the file"
            arg_type = String
            required = true
        "--suffix"
            help = "Suffix to add after chr names, e.g chr1_HS. The '_' is added by default."
            arg_type = String
            default = "HS"
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println(parsed_args)


function get_extension(file)
    return split(file, '.')[end]
end

function modify_gtf(file, suffix)
    output_file = replace(file, r".gtf" => "_sed.gtf")
    open(file, "r") do input
        open(output_file, "w") do output
            for s in eachline(input)
                if !startswith(s, '#')
                    s = split(s, '\t')
                    s[1] = s[1] * "_" * suffix
                    s = join(s, '\t')
                end
                println(output, s)
            end
        end
    end
end

function modify_fasta(file, suffix, extension)
    output_file = replace(file, Regex(".$extension") => "_sed.$extension")
    open(file, "r") do input
        open(output_file, "w") do output
            for s in eachline(input)
                if startswith(s, ">")
                    s = split(s, ' ')
                    s[1] = s[1] * "_" * suffix
                    s = join(s, ' ')
                end
                println(output, s)
            end
        end
    end
end

function main(file, suffix)
    ext = get_extension(file)

    if ext == "gtf"
        modify_gtf(file, suffix)
    elseif ext == "fa" || ext == "fasta"
        modify_fasta(file, suffix, extension)
    else
        error("File extension was not recognized")
    end
end

@time main(parsed_args["file"], parsed_args["suffix"])