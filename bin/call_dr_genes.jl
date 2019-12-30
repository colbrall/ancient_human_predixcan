# call_dr_genes.jl
#
# @author Laura Colbran; Dec. 2019
# functions to call DR status for a group of individuals compared to a reference pop
# assumes models were the same for both populations-- i.e. same genes present in each tissue
#
# julia1.1

using ArgParse
using GZip
using DataFrames
using Statistics

const DELIM = "\t"

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--pop_dir","-d"
            help = "path to directory containing population predictions (1 file per tissue)"
            arg_type = String

        "--sample_list","-s"
            help = "file containing list of samples in pop dir in same order. can be sample.txt file from predixcan run."
            arg_type = String

        "--ref_pop","-r"
            help = "path to directory containing base population to compare to"
            arg_type = String

        "--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String

        "--p_thresh", "-p"
            help = "significance threshold for p-values."
            arg_type = Float64
            default = 0.0
    end
    return parse_args(s)
end

# reads population expression values into a dictionary
function readDict(pop_path::String)
    pop_dict = Dict{String,Array{Float64,1}}()
    GZip.open(pop_path) do f
        e = Array{String,1}()
        for line in eachline(f)
            if startswith(line,"ENSG")
                e = [parse(Float64,d) for d=split(chomp(line),"\t")[2:end]]
                pop_dict[split(chomp(line),"\t")[1]] = e
            end
        end
    end
    return pop_dict
end

# returns array of sample IDs in order from file
function readSamples(path::String)
    s = []
    open(path) do f
        for line in eachline(f)
            append!(s,[split(chomp(line),DELIM)[1]])
        end
    end
    return s
end

# returns tissue name from prediction file
function getTissue(file::String)::String
    return join(split(file,"")[1:end-24])
end

# calculates an empirical 2-sided p-value based on ref population
function popP(stat::Float64,ref::Array{Float64,1})
    #read into dictionary {gene -> pop}
    p = 0.0::Float64
    m = median(ref)
    for pers in ref
        if abs(m - pers) >= abs(m - stat)
            p += 1
        end
    end
    p = p/(length(ref))
    return p,m
end

function drStatus(stat::Float64,ref::Array{Float64,1},thresh::Float64)
    p,m = popP(stat,ref)
    if p <= thresh
        if stat >= m
            return 1
        end
        return -1
    end
    return 0
end

function callDR(pred_dir::String,s_list::String,ref_path::String,out_dir::String,thresh::Float64)
    samples = readSamples(s_list)
    for s in 1:length(samples)
        open("$(out_dir)$(samples[s])_dr_genes.txt","w") do f
            write(f,"#tissue\tgene_id\tdirection\n")
        end
    end
    for item in readdir(pred_dir)
        if !(endswith(item,".full.gz")) continue end
        println("calling $(getTissue(item))...")
        targets = readDict("$(realpath(pred_dir))/$item")
        ref = readDict("$(realpath(ref_path))/$item")
        for gene in keys(targets)
            for s in 1:length(samples)
                dr = drStatus(targets[gene][s],ref[gene],thresh)
                if dr != 0
                    open("$(out_dir)$(samples[s])_dr_genes.txt","a") do f
                        write(f, "$(getTissue(item))\t$(gene)\t$(dr)\n")
                    end
                end
            end
        end
    end
end

function main()
    parsed_args = parseCommandLine()
    callDR(parsed_args["pop_dir"],parsed_args["sample_list"],parsed_args["ref_pop"],parsed_args["out_dir"],parsed_args["p_thresh"])
end

main()
