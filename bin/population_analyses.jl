# comp_pop.jl
# @author Laura Colbran
# functions to make large-scale, population-wide comparisons
#
# contains functions for comparing sets of expression predictions
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
using StatsPlots,Plots,Seaborn
using HypothesisTests

default(color=:blue,leg=false,grid=false,fontfamily="arial",alpha=0.5)
Seaborn.set(style="white")
set_style(Dict("font.family" =>["DejaVu Sans"]))

const DENDRO_SOURCE =
        "/dors/capra_lab/projects/neanderthal_predixcan/bin/WriteDendrogram.R"

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
        "--out"
            help = "to save gene ids passing thresh"
            action = :store_true
        "--thresh","-t"
            nargs = '*'
            help = "correlation(s) above which we'll keep genes. write output file per thresh"
            arg_type = Float64
            default = [0.5]
        "--match_corr"
            action=:store_true
            help = "if you plan on correlating a metric with agreement"

        "--anno","-a"
            arg_type=String
            help="file with individual id and metric you want to correlate to (e.g. output of calc_miss_prop.jl)"
        "--id_map","-m"
            arg_type=String
            help = "sample.txt file with the original ancient ids in the right order for the matched samples."
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
    "uterus" => "uterus", "thyroid" => "thyroid","left_ventricle" => "left_ventricle",
    "Adipose_Subcutaneous" => "Adipose_Subcutaneous","Adipose_Visceral_Omentum" => "Adipose_Visceral_Omentum",
    "Adrenal_Gland" => "Adrenal_Gland","Artery_Aorta" => "Artery_Aorta","Artery_Coronary" => "Artery_Coronary",
    "Artery_Tibial" => "Artery_Tibial","Brain_Amygdala" => "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex" => "Brain_Anterior_cingulate_cortex",
    "Brain_Caudate_basal_ganglia" => "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere" => "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum" => "Brain_Cerebellum","Brain_Cortex" => "Brain_Cortex",
    "Brain_Frontal_Cortex" => "Brain_Frontal_Cortex","Brain_Hippocampus" => "Brain_Hippocampus",
    "Brain_Hypothalamus" => "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia" => "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia" => "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1" => "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra" => "Brain_Substantia_nigra","Breast_Mammary_Tissue" => "Breast_Mammary_Tissue",
    "Cells_Cultured_fibroblasts" => "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes" => "Cells_EBV-transformed_lymphocytes",
    "Colon_Sigmoid" => "Colon_Sigmoid","Colon_Transverse" => "Colon_Transverse",
    "Esophagus_Gastroesophageal_Junction" => "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa" => "Esophagus_Mucosa","Esophagus_Muscularis" => "Esophagus_Muscularis",
    "Heart_Atrial_Appendage" => "Heart_Atrial_Appendage","Heart_Left_Ventricle" => "Heart_Left_Ventricle",
    "Kidney_Cortex" => "Kidney_Cortex",
    "Minor_Salivary_Gland" => "Minor_Salivary_Gland","Muscle_Skeletal" => "Muscle_Skeletal",
    "Nerve_Tibial" => "Nerve_Tibial",
    "Skin_Not_Sun_Exposed_Suprapubic" => "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg" => "Skin_Sun_Exposed_Lower_leg",
    "Small_Intestine_Terminal_Ileum" => "Small_Intestine_Terminal_Ileum",
    "Whole_Blood" => "Whole_Blood","Prostate"=>"Prostate","Testis"=>"Testis","Thyroid"=> "Thyroid",
    "Uterus"=>"Uterus","Vagina"=>"Vagina")
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
    for tiss in keys(map_dict)
        if !in(tiss,["liver","muscle_skeletal","ovary","whole_blood"]) continue end
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
            for id in names(pop)[2:end]
                base_id = split(String(id),"_")[1]
                push!(corr_df,[split(dirname(pop_file),"/")[end],base_id,String(id),tiss,corspearman(temp_pop[id],temp_base[Symbol(base_id)]),cor(temp_pop[id],temp_base[Symbol(base_id)])])
            end
        end
    end
#    comp_plot = Plots.boxplot(corr_df[:comp_pop], corr_df[:rho],notch = true, ylabel="Rho",xlabel= "Comparison Populations",margin=10Plots.mm,ylims=[0,1])
    #tight_layout()
#    Plots.savefig(comp_plot, "pop_corrspearman_box.pdf")

    comp_plot = Seaborn.boxplot(x=corr_df[:comp_pop],y=corr_df[:rho])
    comp_plot = Seaborn.swarmplot(x=corr_df[:comp_pop],y=corr_df[:rho],color=:black,alpha=0.5,size=2)
    comp_plot.set_ylabel("Rho")
    comp_plot.set_xlabel("Comparison Population")
    Seaborn.savefig("pop_corrspearman_box.pdf")
    clf()
    for p in unique(corr_df[:comp_pop])
        println("$p:")
        describe(corr_df[corr_df[:comp_pop] .== p,:rho])
        # comp_plot = Seaborn.boxplot(x=corr_df[:comp_pop],y=corr_df[:rho])
        # comp_plot = Seaborn.swarmplot(x=corr_df[:comp_pop],y=corr_df[:rho],color=:black,alpha=0.5,size=2)
        # comp_plot.set_ylabel("Pred. Reg")
        # comp_plor.set_xlabel("Pred. Reg. $p")
        # Seaborn.savefig("pop_corrspearman_$(p).pdf")
        # clf()
    end
#    comp_plot = Plots.boxplot(corr_df[:comp_pop], corr_df[:r],notch = true, ylabel="Pearson r",xlabel= "Comparison Populations",margin=10Plots.mm)#,ylims=[0,1])
#    Plots.savefig(comp_plot, "pop_corrpearson_box.pdf")

end

#correlation across genes by threshold distribution
function geneSim(base_pop::String,pop_dir::Array{String,1},o::Bool,thresh::Array{Float64,1})
    map_dict = mapNames()::Dict{String,String}
    out = DataFrame(gene=String[],r=Float64[],tissue=String[])
    for tiss in keys(map_dict)
        base = readGeneDF("$(realpath(base_pop))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
        if base == "NA" continue end
        for pop_file in pop_dir
            pop = readGeneDF("$(realpath(pop_file))/$(map_dict[tiss])_elasticNet0_0.5.full.gz")
            if pop == "NA"
                continue
            end
            println("$pop_file, $tiss")
            for i in 1:nrow(base)
                tmp = DataFrame(full=Float64[],sim=Float64[])
                for id in names(pop)[2:end]
                    base_id = split(String(id),"_")[1]
                    if length(pop[pop[:,1] .== base[i,1],id]) == 0 break end
                    tmp = vcat(tmp,DataFrame(full=[base[i,Symbol(base_id)]],sim=pop[pop[:,1] .== base[i,1],id]))
                end
                if nrow(tmp) == 0 continue end
                out = vcat(out,DataFrame(gene=base[i,:gene_id],r=cor(tmp[:full],tmp[:sim]),tissue="$(map_dict[tiss])"))
            end
        end
        #maybe add plot by tiss?
    end
    if o
        for t in thresh
            CSV.write("filtered_genes_$t.txt",out[(out[:r] .> t),:];delim='\t')
        end
    end
    dist_plot = histogram(out[:r],xlabel="r",ylabel="Num. Models",bins=100,margin=10Plots.mm)
    Plots.savefig(dist_plot,"all_gene_r.pdf")
end

function matchCorr(base_pop::String,pop_dir::Array{String,1},miss_prop::String,id_map::String)
    map_dict = mapNames()::Dict{String,String}
    corr_df = DataFrame(comp_pop = String[], base_ind = String[], sim_ind = String[], orig_id = String[],tiss = String[], rho = Float64[],r = Float64[])
    map = CSV.read(id_map;delim='\t',header=[:ind,:fam],allowmissing=:none)[:ind]
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
            i=1
            for id in names(pop)[2:end]
                base_id = split(String(id),"_")[1]
                push!(corr_df,[split(dirname(pop_file),"/")[end],base_id,String(id),map[i],tiss,corspearman(temp_pop[id],temp_base[Symbol(base_id)]),cor(temp_pop[id],temp_base[Symbol(base_id)])])
                i+=1
            end
        end
    end
    println(first(corr_df,3))
    corr_df = join(corr_df,CSV.read(miss_prop;delim='\t',allowmissing=:none),on=:orig_id=>:ind_id)
    corr_df[:all] = "All"
    rho = corspearman(corr_df[:prop_miss],corr_df[:rho])
    println("Spearman rho : $rho; $(OneSampleZTest(atanh(rho), 1, nrow(corr_df)))")
    println("\nCorr summary:")
    describe(corr_df[:rho])
    corr_plot = Seaborn.boxplot(x=corr_df[:tiss],y=corr_df[:rho])
    corr_plot.set_xlabel("Tissue")
    corr_plot.set_ylabel("Rho by Individual")
    Seaborn.savefig("reich_sim_rho_dist_tiss.pdf")
    colours = [:firebrick,:dodgerblue,:darkgreen,:black]
    c=1
    clf()
    for t in unique(corr_df[:tiss])
        corr_plot = Seaborn.kdeplot(corr_df[corr_df[:tiss].== t,:rho],color=colours[c])
        c +=1
    end
    corr_plot.set_xlabel("Rho by Individual")
    Seaborn.savefig("reich_sim_rho_dist_all.pdf")
    clf()
    c=1
    for t in unique(corr_df[:tiss])
        println(t)
        corr_plot = Seaborn.scatter(x=corr_df[corr_df[:tiss].== t,:prop_miss],y=corr_df[corr_df[:tiss].== t,:rho],color=colours[c],alpha=0.25)
        c +=1
    end
    # corr_plot.set_xlabel("Proportion Missing from Genome")
    # corr_plot.set_ylabel("Rho by Individual")
    Seaborn.savefig("reich_sim_rho_miss_corr.pdf")
    clf()
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
        geneSim(parsed_args["base_pop"],parsed_args["pop_dir"],parsed_args["out"],parsed_args["thresh"])
    end
    if parsed_args["match_corr"]
        matchCorr(parsed_args["base_pop"],parsed_args["pop_dir"],parsed_args["anno"],parsed_args["id_map"])
    end
end

main()
