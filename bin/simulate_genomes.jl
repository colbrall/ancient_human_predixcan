# functions to simulate genomes with missing information by several methods
# @author Laura Colbran
#
# julia1.1

using DataFrames
using CSV
using ArgParse
using GZip
using StatsBase

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--random"
            help = "mask random SNPs. Requires -t."
            action = :store_true

        "--matched"
            help = "mask SNPs in pattern matching existing individuals. Requires -m."
            action = :store_true

        "--base_pop","-b"
            help = "path to directory with dosage files for population to start from"
            arg_type = String
            required = true

        "--ids","-i"
            help = "list of IDs from base pop to work from"
            arg_type = String
            required = true

        "--out_dir","-o"
            help = "directory to write output files to"
            arg_type = String
            required = true

        "--num_inds","-n"
            help = "number of genotypes to simulate per threshold or match. Default 1."
            arg_type = Int64
            default = 1

        "--threshold","-t"
            help = "proportion of SNPs to mask, ie 0.1. Default 0.0. Can take multiple."
            nargs = '*'
            arg_type = Float64
            default = [0.0]

        "--matches","-m"
            help = "directory of dosage files for individuals to match to"
            arg_type = String

    end
    return parse_args(s)
end

# updated, expanded version of Maya's simulation
# randomly subsets genotype data to create missingness at a certain threshold.
function matchRandom(base_dir::String,num::Int64,id_path::String,out_dir::String,thresholds::Array{Float64,1})
    ids = CSV.read(id_path,delim='\t';header=["ind_id","index"],allowmissing=:none) #id, and column number of its genotype
    for thresh in thresholds
        picked = ids[sample(axes(ids,1),num;replace=false),:] #choose individuals to downsample
        for i in 1:nrow(picked) #make new ids for simulated individuals
            picked[i,:ind_id] = "$(picked[i,:ind_id])_$thresh"
        end
        run(`mkdir $(realpath(out_dir))/pop_$(thresh)`)
        CSV.write("$(realpath(out_dir))/pop_$(thresh)/sample_$(thresh).txt",DataFrame(x= picked[:ind_id],y=picked[:ind_id]);delim='\t',writeheader=false)
        for item in readdir(base_dir)
            if !endswith(item,".dos.gz") continue end
            out_f = GZip.open("$(realpath(out_dir))/pop_$(thresh)/$(join(split(item,".")[1:(length(split(item,"."))-2)],".")).$(thresh)_missing.dos.gz","w")
            GZip.open("$base_dir/$item","r") do inf
                total = countlines(inf)-1 #pick unique set of SNPs for each simulated individual, minus header
                seek(inf,0)
                snps_to_mask = [sample(collect(1:total),Int64(round(total*thresh));replace=false) for i in 1:nrow(picked)]
                count = 1
                for line in eachline(inf)
                    if startswith(line, "#") continue end
                    info = join(split(line,"\t")[1:6],"\t")
                    pop = split(chomp(line),"\t")[7:end][picked[:index]] # subset pop to picked people
                    # either write dosage as-is for, or swap to '.' if in the line number for them.
                    for i in 1:length(snps_to_mask)
                        if in(count,snps_to_mask[i])
                            pop[i] = "NA"
                        end
                    end
                    write(out_f,"$(info)\t$(join(pop,"\t"))\n")
                    count += 1
                end
            end
        end
    end
end

# subsets genotype data to match pattern of missing SNPs in another set of individuals
# works as expected only if the two have the same SNPs in the same order.
function matchMissing(base_dir::String,id_path::String,out_dir::String,match_dir::String)
    ids = CSV.read(id_path,delim='\t';header=["ind_id","index"],allowmissing=:none) #id, and column number of its genotype
    sampled = false
    picked = DataFrame(ind_id=String[],index=Int64[])
    for item in readdir(base_dir)
        if !endswith(item,".dos.gz") continue end
        chr = split(item,".")[1]
        match_file = filter(f-> startswith(f,chr),readdir(match_dir))[1]
        GZip.open("$base_dir/$item","r") do inf
            GZip.open("$match_dir/$match_file","r") do matchf
                if !sampled #choose individual to downsample
                    num_to_match = length(split(readline(matchf),'\t'))-6
                    picked = ids[sample(axes(ids,1),1),:]
                    picked[1,:ind_id] = "$(picked[1,:ind_id])_match1"
                    for i in 2:num_to_match #make new ids for simulated individuals
                        push!(picked,["$(picked[1,:ind_id])_match$i",picked[1,:index]])
                    end
                    CSV.write("$(realpath(out_dir))/sample_matched.txt",DataFrame(x= picked[:ind_id],y=picked[:ind_id]);delim='\t',writeheader=false)
                    sampled = true
                    seek(matchf,0)
                end
                out_f = GZip.open("$(realpath(out_dir))/$(join(split(item,".")[1:(length(split(item,"."))-2)],".")).matched.dos.gz","w")
                for in_line in eachline(inf)
                    if startswith(in_line, "#") continue end
                    info = join(split(in_line,"\t")[1:6],"\t")
                    #println(info)
                    pop = split(chomp(in_line),"\t")[7:end][picked[:index]]
                    #println(pop)
                    match_line = split(chomp(readline(matchf)),"\t")[7:end]
                    miss_inds = findall(i->(!isequal(i,"0") && !isequal(i,"1") && !isequal(i,"2")),match_line)
                    pop[miss_inds] .= "NA"
                    write(out_f,"$(info)\t$(join(pop,"\t"))\n")
                end
            end
        end
    end
end

function main()
    args = parse_commandline()
    if args["random"]
        println(args["threshold"])
        matchRandom(args["base_pop"],args["num_inds"],args["ids"],args["out_dir"],args["threshold"])
    end
    if args["matched"]
        matchMissing(args["base_pop"],args["ids"],args["out_dir"],args["matches"])
    end
end

main()
