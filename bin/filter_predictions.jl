# filter_predictions.jl
# @author Laura Colbran
# given output of predixcan, filters on gene and/or individual ids
# julia1.1

using DataFrames
using CSV
using ArgParse
using GZip

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--pop_dir","-d"
            help = "path to directory containing population predictions (1 file per tissue)"
            arg_type = String
        "--ind_file","-i"
            help = "list of ind IDs to keep"
            arg_type = String
        "--gene_file","-g"
            help = "list of gene/tissue pairs to keep (output from pop_analyses.jl)"
            arg_type = String
        "--out_dir","-o"
            help = "dir to write filtered output to"
            arg_type = String
    end
    return parse_args(s)
end

function readInds(ind_file::String) #return array of ind IDs as symbols
    inds = Array(Symbol,1)[]
    open(ind_file) do f
        for l in eachline(f)
            append!(inds,[Symbol(split(chomp(l),"\t")[1])])
        end
    end
    return inds
end

function readGenes(gene_file::String,tiss::SubString{String}) #return array of gene IDs to keep in tissue
    genes = CSV.read(gene_file;delim="\t")
    return genes[genes[:tissue] .== tiss,:gene]
end

function filterPreds(pop_dir::String,ind_file,gene_file,out_dir::String)
    for file in readdir(pop_dir)
        if !endswith(file, "full.gz") continue end
        tiss = split(file,"_elasticNet0")[1]
        println(tiss)
        tmp = CSV.read(GZip.open("$pop_dir$file");delim='\t',allowmissing=:none)
        if typeof(ind_file) != Nothing
            tmp = tmp[:,readInds(ind_file)]
        end
        if typeof(gene_file) != Nothing
            targets = readGenes(gene_file,tiss)
            tmp = tmp[[in(i,targets) for i in tmp[:gene]],:]
        end
        outf = "$(out_dir)$(splitext(file)[1])"
        CSV.write(outf,tmp;delim="\t")
        run(`gzip $outf`)
    end
end

function main()
    parsed_args = parse_commandline()
    filterPreds(parsed_args["pop_dir"],parsed_args["ind_file"],parsed_args["gene_file"],parsed_args["out_dir"])
end

main()
