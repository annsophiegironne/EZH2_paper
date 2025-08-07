using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--dir"
            help = "Where the files to merge are stored"
            arg_type = String
            required = true
        "--skip_dir"
            help = "Directory(ies) to filter out when fetching count files"
            arg_type = String
            nargs = '*'
        "--counts"
            help = "Column where to find counts to merge"
            required = true
            arg_type = UInt8
            dest_name = "counts_column"
        "--sep"
            help = "Separator of the input files"
            arg_type = String
            default = "\t"
        "--gene_names"
            help = "Whether to add gene names or not"
            arg_type = Bool
            default = true
            dest_name = "add_gene_names"
        "--pattern"
            help = "Pattern to get all the files to merge"
            arg_type = String
            default = "ReadsPerGene.out.tab"
        "--skip_pattern"
            help = "Pattern to filter out genes. Must have \"\", \"ENST|ENSMUST\" for example."
            arg_type = Regex
            default = "nothing"
            dest_name = "pattern_to_filter_out"
        "--delim"
            help = "Delimiter of the output file"
            default = ","
            arg_type = String
        "--ids"
            help = "Column where to find gene ids"
            default = 1
            arg_type = UInt8
            dest_name = "ids_column"
        "--header"
            help = "Length of header"
            default = 4
            arg_type = UInt8
            dest_name = "length_of_header"
        "--output"
            help = "Name of the output file. If directory not specified in name, output file will be written in current directory"
            default = "./merged_counts.csv"
            dest_name = "output_file"
            arg_type = String
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println(parsed_args)

function list_files(dir, pattern, skip_dir)
    files_list = Vector()
    for (root, dirs, files) in walkdir(dir)
        for file in files
            if contains(file, pattern)
                file = root * "/" * file
                push!(files_list, file)
            end
        end
    end

    files_list = filter(x -> !any(occursin.(skip_dir, x)), files_list)
    return files_list
end

function get_ids(files_list, sep, ids_column, length_of_header, pattern_to_filter_out)
    file::IOStream = open(files_list[1], "r")::IOStream
    ids = Vector{String}()
    line = 1
    while !eof(file)
        if line <= length_of_header
            s = readline(file)
            line += 1
        else
            s = readline(file)
            s = split(s, sep)[ids_column]
            if !occursin(pattern_to_filter_out, s)
                push!(ids, s)
            end
            line += 1
        end
    end
    return ids
end

function merge_files(sep, ids_column, counts_column, files_to_merge, ids, length_of_header, pattern_to_filter_out)
    m = length(ids)
    n = length(files_to_merge)
    counts = Array{Any}(missing, m, 1)
    counts[:,1] = @view(ids[:])
    sampleNames::Vector{String} = Vector{String}(undef, n)

    for i in 1:length(files_to_merge)
        file::IOStream = open(files_to_merge[i], "r")::IOStream
        line = 1
        lcounts = Vector{String}(undef, m)
        filtered_counts = 0
        while !eof(file)
            if line <= length_of_header
                s = readline(file)
                line += 1
            else
                s = readline(file)
                s = split(s, sep)
                if !occursin(pattern_to_filter_out, s[ids_column])
                    filtered_counts += 1
                    if s[ids_column] == ids[filtered_counts]
                        lcounts[filtered_counts] = s[counts_column]
                    else
                        error("Files need to be in the same order")
                    end
                end
                line += 1
            end
        end
    sampleName = split(files_to_merge[i], "/")[end-1]
    sampleNames[i] = sampleName
    counts = hcat(counts, lcounts)
    close(file)
    end
    final_matrix = vcat(permutedims(vcat(["ids"], sampleNames)), counts)
    return final_matrix
end

function add_gene_names_to_matrix(matrix, gencode)
    file = open("annotations/Homo_sapiens/Homo_sapiens.GRCh38.Gencode40.mapping.genes.tsv", "r")
    dict = Dict{String, String}()
    while !eof(file)
        s = readline(file)
        s = split(s, "\t")
        dict[s[1]] = s[2]
    end
    rows = size(matrix)[1]
    add = Vector{Any}(undef, rows)
    for i in 1:rows
        if i == 1
            add[i] = "gene_id"
        else
            if haskey(dict, matrix[i, 1])
                add[i] = dict[matrix[i, 1]]
            else
                add[i] = matrix[i, 1]
            end
        end
    end
    matrix = hcat(matrix, add)
    return(matrix)
end

function write_matrix_to_file(matrix, output_file, delim)
    file = open(output_file, "w")
    rows = size(matrix)[1]

    for i in 1:rows
        println(file, join(matrix[i,:], delim))
    end

    close(file)
end

function run(pattern, sep, delim, counts_column, ids_column, dir, output_file, length_of_header, pattern_to_filter_out, add_gene_names, gencode, skip_dir)
    files_list = list_files(dir, pattern, skip_dir)
    ids = get_ids(files_list, sep, ids_column, length_of_header, pattern_to_filter_out)
    matrix = merge_files(sep, ids_column, counts_column, files_list, ids, length_of_header, pattern_to_filter_out)
    if add_gene_names == true
        final_matrix = add_gene_names_to_matrix(matrix, gencode)
        write_matrix_to_file(final_matrix, output_file, delim)
    else
        write_matrix_to_file(matrix, output_file, delim)
    end
end

run(parsed_args["pattern"], parsed_args["sep"], parsed_args["delim"], parsed_args["counts_column"], 
    parsed_args["ids_column"], parsed_args["dir"], parsed_args["output_file"], parsed_args["length_of_header"],
    parsed_args["pattern_to_filter_out"], parsed_args["add_gene_names"], parsed_args["skip_dir"])
