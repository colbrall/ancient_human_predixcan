# comp_pop.jl
# @author Laura Colbran
# functions to compare PrediXcan results in archaics to Eric's bioVU results or 1kG
#
# contains functions for calculating empirical p-values, calling DR genes, direction bias, PCA,
# and computing a distance matrix for a tree
#
# runs with Julia version 1.1 and requires R3.3 or greater to be present in the environment with package dendextend

using ArgParse
using GZip
using DataFrames
using CSV
using Distances
using RCall
using StatsBase
using Statistics
using StatsPlots
using Plots

default(color=:blue,leg=false,grid=false,fontfamily="arial",alpha=0.5)
#using PyPlot

const DENDRO_SOURCE =
        "/dors/capra_lab/projects/neanderthal_predixcan/bin/WriteDendrogram.R"
#const N = 18621 #number people in population for empirical p-value calculation
const N = 2504

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--pop_dir","-d"
            nargs= '*'
            help = "path(s) to directory containing population predictions (1 file per tissue)"
            arg_type = String
        "--pop_ids","-i"
            help = "path to file with IDs for all individuals being considered"
            arg_type = String
        "--filter","-f"
            help = "path to gene list to filter on for tree-building. optional."
            arg_type = String
        "--exclude"
            help = "whether to exclude genes (as opposed to include only those). Default to false."
            action = :store_true
        "--spearman"
            help = "whether to use spearman correlation for distance matrix (rather than pearson). Default to false."
            action = :store_true
        "--tree"
            help = "hierarchical clustering to return newick tree. requires --pop_dir --pop_ids and WriteDendrogram.R"
            action = :store_true
        "--base_pop","-b"
            help = "path to directory containing base population to compare to"
            arg_type = String
        "--individuals"
            help = "to compute similarities between individuals"
            action = :store_true
        "--genes"
            help = "to compute similarities between genes"
            action = :store_true
    end
    return parse_args(s)
end

# reads gene x tiss matrix into DataFrame
function readGeneDF(f_path::String)
  if ispath(f_path)
    df = CSV.read(GZip.open(f_path);delim='\t',allowmissing=:none)
    try rename!(df, :gene=> :gene_id) finally return df end
  else
    df = "NA"
    return df
  end
end

# reads file into DataFrame
function readDF(f_path::String)
  if ispath(f_path)
    df = CSV.read(GZip.open(f_path);delim='\t',allowmissing=:none, normalizenames=false)
  else
    df = "NA"
  end
  return df
end

# reads population expression values into a dictionary
function readDict(path::String)
  dict = Dict{String,Array{SubString{String},1}}()
  try
    GZip.open(path) do f
      for line in eachline(f)
        if startswith(line,"ENSG")
          e = split(chomp(line),"\t")
          dict[e[1]] = e[2:end]
        end
      end
    end
    return dict
  catch
    return "NA"
  end
end

# reads population and ids into a dictionary
function popDict(path::String)
  dict = Dict{SubString{String},SubString{String}}()
  try
    GZip.open(path) do f
      for line in eachline(f)
        e = split(chomp(line),"\t")
        dict[e[1]] = e[2]
      end
    end
    return dict
  catch
    return "NA"
  end
end
# my tissue names -> Eric's tissue names
function mapNames()
  return Dict{String,String}(
    "adipose_subcutaneous" => "Adipose-Subcutaneous",
    "adipose_visceral_omentum" => "adipose_visceral_omentum",
    "brain_putamen_basal_ganglia" => "Brain-Putamen-basalganglia",
    "pancreas" => "Pancreas", "breast_mammary_tissue" => "Breast-MammaryTissue",
    "pituitary" => "Pituitary", "adrenal_gland" => "AdrenalGland",
    "cells_ebv_transformed_lymphocytes" => "cells_ebv-transformed_lymphocytes",
    "anterior_cingulate_cortex" => "Brain-Anteriorcingulatecortex-BA24",
    "cells_transformed_fibroblasts" => "Cells-Transformedfibroblasts",
    "skin_nosun_suprapubic" => "Skin-NotSunExposed-Suprapubic",
    "artery_aorta" => "Artery-Aorta", "colon_sigmoid" => "Colon-Sigmoid",
    "skin_sun_lower_leg" => "Skin-SunExposed-Lowerleg",
    "artery_coronary" => "Artery-Coronary", "colon_transverse" => "Colon-Transverse",
    "small_intestine_terminal_ileum" => "SmallIntestine-TerminalIleum",
    "artery_tibial" => "Artery-Tibial", "spleen" => "Spleen",
    "esophagus_gastroesophageal_junction" => "Esophagus-GastroesophagealJunction",
    "brain_caudate_basal_ganglia" => "Brain-Caudate-basalganglia",
    "esophagus_mucosa" => "Esophagus-Mucosa", "stomach" => "Stomach",
    "brain_cerebellar_hemisphere" => "Brain-CerebellarHemisphere",
    "esophagus_muscularis" => "Esophagus-Muscularis",
    "brain_cerebellum" => "Brain-Cerebellum", "nerve_tibial" => "Nerve-Tibial",
    "heart_atrial_appendage" => "Heart-AtrialAppendage",
    "brain_cortex" => "Brain-Cortex", "liver" => "Liver", "lung" => "Lung",
    "brain_frontal_cortex" => "Brain-FrontalCortex-BA9",
    "brain_hippocampus" => "Brain-Hippocampus", "muscle_skeletal" => "Muscle-Skeletal",
    "whole_blood" => "WholeBlood", "brain_hypothalamus" => "Brain-Hypothalamus",
    "brain_nucleus_accumbens_basal_ganglia" => "Brain-Nucleusaccumbens-basalganglia",
    "prostate" => "prostate","ovary" => "Ovary", "vagina" => "vagina", "testis" => "testis",
    "uterus" => "uterus", "thyroid" => "thyroid","left_ventricle" => "left_ventricle")
end

# reverses key:value dictionary; naive-- doesn't check for unique values
function flipDict(dict::Dict{String,String})
  new_dict = Dict{String,String}()
  for key in keys(dict)
    new_dict[dict[key]] = key
  end
  return new_dict
end

#1kG pop -> superpop
function thousGenSuperPopDict()
  return Dict{String,String}(
    "GWD" => "AFR", "MSL" => "AFR","ESN" => "AFR", "MXL" => "AMR", "CLM" => "AMR",
    "PEL" => "AMR", "TSI" => "EUR", "IBS" => "EUR", "PJL" => "SAS", "STU" => "SAS",
    "ITU" => "SAS", "GBR" => "EUR", "CHB" => "EAS", "JPT" => "EAS", "CDX" => "EAS",
    "YRI" => "AFR", "LWK" => "AFR", "CEU" => "EUR", "GIH" => "SAS", "ASW" => "AFR",
    "CHS" => "EAS", "KHV" => "EAS", "ACB" => "AFR", "PUR" => "AMR", "BEB" => "SAS",
    "FIN" => "EUR")
end


# assembles distance matrix
function distMat(pop_dir::Array{String,1},filter,to_excl::Bool,spear::Bool,pop_dict::Dict{SubString{String},SubString{String}},tiss::String)
    excl = ["MXL", "CLM", "PUR", "ACB","ASW", "PJL", "PEL"]
    map_dict = mapNames()::Dict{String,String}
    if length(pop_dir) == 0 #take first population, make it be base DF
      df = "none"
      println("No Population Directory")
    else
      df = readGeneDF("$(realpath(pop_dir[1]))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
    end
    if df == "NA" return df end
    # get values for tissue from other pop directories
    if (length(pop_dir) > 1)
        for p in pop_dir[2:end]
            tmp = readGeneDF("$(realpath(p))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
            if tmp != "NA" df = join(tmp, df, on = :gene_id) end
        end
    end
    # remove people to exclude (admixed populations)
    for ind in keys(pop_dict)
      if in(pop_dict[ind], excl)
          try deletecols!(df,Symbol(ind)) finally continue end
      end
    end
    # filter gene list
    if typeof(filter) != Nothing
        filter_genes = CSV.read(filter, delim='\t',header=["gene_id"])
        if to_excl #ie you want to exclude the genes in your filtered list
            df = join(df,filter_genes,on=:gene_id,kind=:anti)
        else
            df = join(df,filter_genes,on=:gene_id,kind=:inner)
        end
    end
    deletecols!(df, :gene_id)
    if spear
        for col in names(df)
            df[col] = tiedrank(df[col]) #convert to ranks so CorrDist is actually Spearman correlation.
        end
    end
    #writedlm("$(tiss)_distance.txt",pairwise(CorrDist(),Matrix(df)))
    return pairwise(CorrDist(),Matrix(df),dims=2)
end

# hierarchical clustering on a distance matrix, saves newick tree
function makeTree(id_path::String,pop_dir::Array{String,1},filter,to_excl::Bool,spear::Bool)
    println("Clustering.....")
    map_dict = mapNames()::Dict{String,String}
    for key in keys(map_dict)
        pop_dict = popDict(id_path)
        dist = distMat(pop_dir::Array{String,1},filter,to_excl::Bool,spear::Bool,pop_dict,key)
        if dist == "NA" continue end
        out_tree = "$(key)_dendrogram.newick"
        R"""
            suppressPackageStartupMessages(library(dendextend))
            source($DENDRO_SOURCE)
            dend <- as.dendrogram(hclust(as.dist($dist)), hang = -1)
            l <- order.dendrogram(dend)
            WriteDendrogram(dend,file=$out_tree,quoteLabels = F)
        """
    end
end

#correlation across individuals by threshold-- boxplots
function indSim(base_pop::String,pop_dir::Array{String,1})
    map_dict = mapNames()::Dict{String,String}
    corr_df = DataFrame(comp_pop = String[], base_ind = String[], sim_ind = String[], tiss = String[], rho = Float64[],r = Float64[])
    uni_df = DataFrame(id = String[],comp_pop = String[],num_comp_uniq=Int64[],num_base_uniq=Int64()) #id = model_tiss
    for tiss in keys(map_dict)
        base = readGeneDF("$(realpath(base_pop))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
        if base == "NA" continue end
        for pop_file in pop_dir
            pop = readGeneDF("$(realpath(pop_file))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
            if pop == "NA"
                continue
            end
            println("$pop_file, $tiss")
            temp_base = sort!(join(base,pop,on=:gene_id,kind=:semi),cols=[:gene_id])
            temp_pop = sort!(join(pop,base,on=:gene_id,kind=:semi),cols=[:gene_id])
            for i in 1:nrow(temp_pop)
                push!(uni_df,["$(temp_pop[i,:gene_id])_$tiss",split(dirname(pop_file),"/")[end],length(unique(temp_pop[i,2:end])),length(unique(temp_base[i,2:end]))])
            end
            for id in names(pop)[2:end]
                base_id = split(String(id),"_")[1]
                push!(corr_df,[split(dirname(pop_file),"/")[end],base_id,String(id),tiss,corspearman(temp_pop[id],temp_base[Symbol(base_id)]),cor(temp_pop[id],temp_base[Symbol(base_id)])])
            end
        end
    end
    comp_plot = boxplot(corr_df[:comp_pop], corr_df[:rho],notch = true, ylabel="Rho",xlabel= "Comparison Populations",margin=10Plots.mm)
    #tight_layout()
    Plots.savefig(comp_plot, "pop_corrspearman_box.pdf")
    comp_plot = boxplot(corr_df[:comp_pop], corr_df[:r],notch = true, ylabel="Pearson r",xlabel= "Comparison Populations",margin=10Plots.mm)
    Plots.savefig(comp_plot, "pop_corrpearson_box.pdf")
    comp_plot = violin(uni_df[:comp_pop],uni_df[:num_base_uniq],side=:left)
    comp_plot = violin!(uni_df[:comp_pop],uni_df[:num_comp_uniq],side=:right,color=:red)
    Plots.savefig(comp_plot, "num_uniq_values.pdf")
end

#correlation across genes by threshold distribution
function geneSim(base_pop::String,pop_dir::Array{String,1})
    map_dict = mapNames()::Dict{String,String}
    out = DataFrame(gene=String[],rho=Float64[],tissue=String[])
    for tiss in keys(map_dict)
        println("$tiss")
        base = readGeneDF("$(realpath(base_pop))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
        for pop_file in pop_dir
            pop = readGeneDF("$(realpath(pop_file))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
            for i in 1:nrow(base)
                tmp = DataFrame(full=Float64[],sim=Float64[])
                for id in names(pop)[2:end]
                    base_id = split(String(id),"_")[1]
                    tmp = vcat(tmp,DataFrame(full=[base[i,Symbol(base_id)]],sim=[pop[i,id]]))
                end
                out = vcat(out,DataFrame(gene=base[i,:gene_id],rho=cor(tmp[:full],tmp[:sim]),tissue=tiss))
            end
        end
        #maybe add plot by tiss?
    end
    println(first(out,10))
    println(nrow(out))
    dist_plot = histogram(out[:rho],xlabel="Rho",ylabel="Num. Models",bins=100,margin=10Plots.mm)
    Plots.savefig(dist_plot,"all_gene_rho.pdf")
end

function main()
    parsed_args = parse_commandline()
#### Clustering
    if parsed_args["tree"]
        makeTree(parsed_args["pop_ids"],parsed_args["pop_dir"],parsed_args["filter"],parsed_args["exclude"],parsed_args["spearman"])
    end
    if parsed_args["individuals"]
        indSim(parsed_args["base_pop"],parsed_args["pop_dir"])
    end
    if parsed_args["genes"]
        geneSim(parsed_args["base_pop"],parsed_args["pop_dir"])
    end
end

main()
