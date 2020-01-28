# compare_by_group.jl
# @author Laura Colbran
#
# compare predictions between groups of individuals
#
# julia1.1

using ArgParse
using GZip
using DataFrames
using CSV
using Statistics
using StatsBase
using HypothesisTests
using Plots
using Seaborn

default(color=:black,leg=false,grid=false,fontfamily="arial",alpha=0.5)
Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--pop_dir","-d"
            help = "path to directory containing population predictions (1 file per tissue). required. can be file path for plotting specific prediction"
            arg_type = String

        "--sample_list","-s"
            help = "file containing sample IDs and group annotations. required."
            arg_type = String

        "--column","-c"
            help = "column of sample file that contains groups/variable you want to compare across. required."
            arg_type = Int64

        "--groups","-g"
            nargs='*'
            help = "groups to compare. must match designation in sample_list. if blank, assumes you want continuous comparison across the column"
            arg_type = String

        "--multi_correct","-m"
            help = "type of multiple testing correction to do across genes. can be none, bonferroni, fdr"
            arg_type = String
            default = "none"

        "--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"

        "--targets","-t"
            help = "bed file with specific genes to check. if empty, does all genes."
            arg_type = String

        "--plot"
            help = ""
            action=:store_true

        "--genes","-n"
            nargs='*'
            help = "gene id(s) to plot"
            arg_type = String
    end
    return parse_args(s)
end

# returns bonferroni multiple testing correction significance threshold
function bonferroni(m::Int64)
  return 0.05/m
end

# returns Benjamini Hochberg FDR multiple testing correction significance threshold
function FDR(a::Array{Float64,1})
  fdr = 0.0
  f = DataFrame(values = sort(a), rank = collect(1:length(a)))
  for i in 1:nrow(f)
    if (0.05*f[i,2]/nrow(f)) < f[i,1]
      fdr = f[i,1]
      break
    end
  end
  return fdr
end

# returns dict of group=>samp_ids
function sampDict(s_path::String,groups::Array{String,1},col::Int64)
    dict = Dict{SubString,Array{SubString,1}}()
    for g in groups
        dict[g] = Array{SubString,1}()
    end
    dict["others"] = Array{SubString,1}()
    open(s_path) do f
        for line in eachline(f)
            l = split(chomp(line),'\t')
            if startswith(line,"#")
                println("$(l[col]) groups to compare:")
                continue
            end
            try
                append!(dict[l[col]],[l[1]])
            catch
                append!(dict["others"],[l[1]])
            end
        end
    end
    if length(groups) != 1
        delete!(dict,"others")
    end
    for g in keys(dict)
        println("$g: $(length(dict[g]))")
    end
    println("")
    return dict
end

function sampDF(s_path::String,col::Int64)
    df = CSV.read(s_path; delim='\t',allowmissing=:none)
    return df[:,[1,col]]
end

function pullInds(d::Dict{SubString,Array{SubString,1}},l::Array{SubString{String},1})
    dict = Dict{SubString,Array{Int64,1}}()
    for g in keys(d)
        dict[g] = [findfirst(x->x==i,l) for i in d[g]]
    end
    return dict
end

# returns tissue name from prediction file
function getTissue(file::String)::String
    return join(split(file,"")[1:end-24])
end

function geneTargets(file)::Array{SubString,1}
    arr = SubString[]
    if typeof(file) == Nothing return arr end
    open(file) do f
        for line in eachline(f)
            if startswith(line,"#") continue end
            append!(arr,[split(chomp(line),"\t")[4]])
        end
    end
    return arr
end

function kwTest(arr::Array{Array{Float64,1},1})
    p = 1.0
    if length(arr) == 3
        p = pvalue(KruskalWallisTest(arr[1],arr[2],arr[3]))
    elseif length(arr) == 4
        p = pvalue(KruskalWallisTest(arr[1],arr[2],arr[3],arr[4]))
    elseif length(arr) == 5
        p = pvalue(KruskalWallisTest(arr[1],arr[2],arr[3],arr[4],arr[5]))
    else
        println("ERROR: are you sure you want to compare more than 5 groups? add an elseif to kwTest() if so.")
        exit()
    end
    return p
end

function corrVar(pop_path::String,samp_path::String,col::Int64,correct::String,out_path::String,targ_file::String)
    var = sampDF(samp_path,col)
    gene_res = DataFrames.DataFrame(gene=String[],tiss=String[],rho=Float64[],p=Float64[])
    targets = geneTargets(targ_file)
    for file in readdir(pop_path)
        if !endswith(file,"full.gz") continue end
        println(file)
        indices = Dict{SubString,Array{Int64,1}}()
        GZip.open("$(pop_path)$file") do f
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene")
                    indices = [findfirst(x->x==i,l) for i in var[:,1]]
                    continue
                end
                if length(targets) > 0 && !in(l[1],targets) continue end
                rho = corspearman(var[:,2],parse.(Float64,l[indices]))
                push!(gene_res,[l[1],getTissue(file),rho,pvalue(OneSampleZTest(atanh(rho), 1, nrow(var)))])
            end
        end
        thresh = 0.0
        if correct == "bonferroni"
            thresh = bonferroni(nrow(gene_res))
            gene_res[:bonf_pass] = false
        elseif correct == "fdr"
            thresh = FDR(gene_res[:p])
            gene_res[:fdr_pass] = false
        else
            thresh = 0.05
            gene_res[:uncorr_pass] = false
        end
        for row in 1:nrow(gene_res)
            if gene_res[row,:p] < thresh
                gene_res[row,ncol(gene_res)] = true
            end
        end
        println("genes passing $correct correction: $(nrow(gene_res[gene_res[:,ncol(gene_res)] .== true,:]))/$(nrow(gene_res))")
        println("threshold = $thresh\n")
    end
    if isfile("$(realpath(out_path))/pvalues_$(join(group_names,"_")).txt")
        CSV.write("$(realpath(out_path))/pvalues_$(join(group_names,"_")).txt",gene_res;delim='\t',append=true)
    else
        CSV.write("$(realpath(out_path))/pvalues_$(join(group_names,"_")).txt",gene_res;delim='\t')
    end
end

function grpComp(pop_path::String,samp_path::String,group_names::Array{String,1},col::Int64,correct::String,out_path::String,targ_file)
    groups = sampDict(samp_path,group_names,col)
    targets = geneTargets(targ_file)
    gene_res = DataFrames.DataFrame()
    for file in readdir(pop_path)
        if !endswith(file,"full.gz") continue end
        println(file)
        indices = Dict{SubString,Array{Int64,1}}()
        GZip.open("$(pop_path)$file") do f
            gene_res = DataFrames.DataFrame(gene=String[],tiss=String[],p=Float64[])
            for g in keys(groups)
                gene_res[Symbol("$(g)_med")] = Float64[]
            end
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene")
                    indices = pullInds(groups,l)
                    continue
                end
                if length(targets) > 0 && !in(l[1],targets) continue end
                newrow = [l[1],getTissue(file),-1]
                append!(newrow,fill(-1,length(keys(groups))))
                push!(gene_res, newrow)
                comp = Array{Float64,1}[]
                for g in keys(groups)
                    append!(comp,[parse.(Float64,l[indices[g]])])
                    gene_res[nrow(gene_res),Symbol("$(g)_med")] = median(parse.(Float64,l[indices[g]])) #should be 0 if same as gtex
                end
                if length(keys(groups)) == 2
                    gene_res[nrow(gene_res),:p] = pvalue(MannWhitneyUTest(comp[1],comp[2]))
                else
                    gene_res[nrow(gene_res),:p] = kwTest(comp)
                end
            end
        end
        thresh = 0.0
        if correct == "bonferroni"
            thresh = bonferroni(nrow(gene_res))
            gene_res[:bonf_pass] = false
        elseif correct == "fdr"
            thresh = FDR(gene_res[:p])
            gene_res[:fdr_pass] = false
        else
            thresh = 0.05
            gene_res[:uncorr_pass] = false
        end
        for row in 1:nrow(gene_res)
            if gene_res[row,:p] < thresh
                gene_res[row,ncol(gene_res)] = true
            end
        end
        println("genes passing $correct correction: $(nrow(gene_res[gene_res[:,ncol(gene_res)] .== true,:]))/$(nrow(gene_res))")
        println("threshold = $thresh\n")
        if isfile("$(realpath(out_path))/pvalues_$(join(group_names,"_")).txt")
            CSV.write("$(realpath(out_path))/pvalues_$(join(group_names,"_")).txt",gene_res;delim='\t',append=true)
        else
            CSV.write("$(realpath(out_path))/pvalues_$(join(group_names,"_")).txt",gene_res;delim='\t')
        end
    end
end

function plotSwarm(pop_path::String,samp_path::String,group_names::Array{String,1},col::Int64,out_path::String,gene_ids::Array{String,1})
    groups = sampDict(samp_path,group_names,col)
    gene_res = DataFrames.DataFrame()
    if isfile(pop_path)
        GZip.open("$(pop_path)") do f
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene")
                    gene_res[:samp_id] = l[2:end]
                    gene_res[:groups] = "Other"
                    indices = pullInds(groups,l)
                    # println(l[indices["Agri."]])
                    # println(gene_res[indices["Agri."].-1,:samp_id])
                    for g in keys(groups)
                        gene_res[indices[g].-1,:groups] = g
                    end
                    continue
                end
                if !in(l[1], gene_ids) continue end
                gene_res[Symbol(l[1])] = [parse(Float64,i) for i in l[2:end]]
                #CSV.write("$(l[1])_$(join(group_names,'_')).txt",gene_res[[in(i,group_names) for i in gene_res[:groups]],:];delim='\t')
                s_plot = swarmplot(gene_res[[in(i,group_names) for i in gene_res[:groups]],:groups],
                                gene_res[[in(i,group_names) for i in gene_res[:groups]],Symbol(l[1])])
                s_plot.set_title("$(l[1])")
                s_plot.set_ylabel("Pred. Norm. Expr.")
                Seaborn.savefig("$(l[1])_$(join(group_names,'_')).pdf")
                clf()
            end
        end
    else
        println("not implemented yet-- please specify tissue")
        #plot all tiss
    end
end

function plotTrend(pop_path::String,samp_path::String,col::Int64,out_path::String,gene_ids::Array{String,1})
    var = sampDF(samp_path,col)
    if isfile(pop_path)
        indices = Dict{SubString,Array{Int64,1}}()
        GZip.open("$(pop_path)") do f
            for line in eachline(f)
                l = split(chomp(line),'\t')
                if startswith(line,"gene")
                    indices = [findfirst(x->x==i,l) for i in var[:,1]]
                    continue
                end
                if !in(l[1], gene_ids) continue end
                var[Symbol(l[1])] = parse.(Float64,l[indices])
                s_plot = Plots.scatter(var[:,2],var[Symbol(l[1])],xlabel="$(names(var)[2])",ylabel="Pred. Expr.",title ="$(l[1])",margin=10Plots.mm)
                Plots.savefig(s_plot,"$(l[1])_$(names(var)[2]).pdf")
            end
        end
    else
        println("not implemented yet-- please specify tissue")
        #plot all tiss
    end
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["plot"] && length(parsed_args["groups"]) != 0
        plotSwarm(parsed_args["pop_dir"],parsed_args["sample_list"],parsed_args["groups"],parsed_args["column"],(parsed_args["out_dir"]),parsed_args["genes"])
    elseif parsed_args["plot"] && length(parsed_args["groups"]) == 0
        plotTrend(parsed_args["pop_dir"],parsed_args["sample_list"],parsed_args["column"],(parsed_args["out_dir"]),parsed_args["genes"])
    elseif length(parsed_args["groups"]) == 0
        corrVar(parsed_args["pop_dir"],parsed_args["sample_list"],parsed_args["column"],parsed_args["multi_correct"],(parsed_args["out_dir"]),parsed_args["targets"])
    else #if 1, compares that group to everyone else
        grpComp(parsed_args["pop_dir"],parsed_args["sample_list"],parsed_args["groups"],parsed_args["column"],parsed_args["multi_correct"],(parsed_args["out_dir"]),parsed_args["targets"])
    end
end

main()
