# @author Laura Colbran
# 2019-02-20
# compares how well two sets of PrediXcan models match
#
# julia 1.1

using ArgParse
using DataFrames
using StatsBase
using HypothesisTests
using StatsPlots
using SQLite
using GZip

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--databases","-d"
            nargs='*'
            help = "paths to the model databases"
            arg_type = String
            required = true
        "--performance"
            help = "to run comparison between certain models' performances and SNP counts across different training instances"
            action = :store_true
        # "--snps"
        #     help = "to compare specific SNPs and their weights"
        #     action = :store_true
        "--genomes","-g"
            help = "path to directory containing genomes"
            arg_type = String
        "--missingness"
            help = "to do analyses on missing SNPs by model"
            action = :store_true
    end
    return parse_args(s)
end

#compares R2 and Number of SNPs in models from both sets
function performance(db_files::Array{String,1})
    q = "SELECT * FROM 'extra'"
    models_1 = DataFrame(SQLite.Query(SQLite.DB(db_files[1]),q))
#    models_1 = models_1[models_1[Symbol("pred.perf.R2")] .> 0.1,:]
    models_2 = DataFrame(SQLite.Query(SQLite.DB(db_files[2]),q))
    println("Number of models in $(basename(db_files[1])): $(nrow(models_1))")
    println("Number of models in $(basename(db_files[2])): $(nrow(models_2))")
    all_models = join(models_1,models_2,on=:gene,kind=:inner,makeunique=true)
    println("Number models shared: $(nrow(all_models))")
    #println(join(models_2,models_1,on=:gene,kind=:anti,makeunique=true))
    #exit()
    #println(first(all_models,6))
    println("R2:")
        rho = corspearman(convert(Array{Float64,1},all_models[Symbol("pred.perf.R2")]),convert(Array{Float64,1},all_models[Symbol("pred.perf.R2_1")]))
	println("In order of arguments:")
	describe(all_models[Symbol("pred.perf.R2")])
	describe(all_models[Symbol("pred.perf.R2_1")])
        #rho = corspearman(convert(Array{Float64,1},all_models[Symbol("pred.perf.R2")]),convert(Array{Float64,1},all_models[Symbol("R2")]))
        println("Spearman rho: $(rho)")
        println(OneSampleZTest(atanh(rho), 1, nrow(all_models)))
    println("Num SNPs:")
        rho = corspearman(convert(Array{Float64,1},all_models[Symbol("n.snps.in.model")]),convert(Array{Float64,1},all_models[Symbol("n.snps.in.model_1")]))
        #rho = corspearman(convert(Array{Float64,1},all_models[Symbol("n.snps.in.model")]),convert(Array{Float64,1},all_models[Symbol("n.snps")]))
        println("Spearman rho: $(rho)")
        println(OneSampleZTest(atanh(rho), 1, nrow(all_models)))
end

# for each SNP, count the number of people missing it, and for each individual,
function missingDict(path::String)
    ind_dict = Dict{String,Array{String,1}}()
    snp_dict = Dict{String,Array{String,1}}()
    for item in readdir(path)
        if !endswith(item,"dos.gz") continue end
        println(item)
        GZip.open("$(realpath(path))/$item") do f
            for line in eachline(f)
                l = split(chomp(line),"\t")
                ind_ids = [string(i) for i in collect(1:length(l)-6)] #if ids aren't present, ids will just be index of individual from 1
                if startswith(line,"#")
                    ind_ids = l[7:length(l)]
                    continue
                end
                for i in collect(7:length(l))
                    if l[i] == "NA"
                        #println(line)
                        if in(l[2], keys(snp_dict))
                            append!(snp_dict[l[2]],[ind_ids[i-6]])
                        else
                            snp_dict[l[2]] = [ind_ids[i-6]]
                        end
                        if in(ind_ids[i-6], keys(ind_dict))
                            append!(ind_dict[ind_ids[i-6]],[l[2]])
                        else
                            ind_dict[ind_ids[i-6]] = [l[2]]
                        end
                    end
                end
            end
        end
    end
    return ind_dict,snp_dict
end

# analyses of missing SNPs in models
function missingSNPs(db_files::Array{String,1},genome_path::String)
    ind_dict,snp_dict = missingDict(genome_path)
    for db in db_files
        name = splitext(splitpath(db)[end])[1]
        println(name)
        gene_snps = Dict{String,Array{String,1}}()
        q = "SELECT * FROM 'weights'"
        for snp in SQLite.Query(SQLite.DB(db),q)
            #println(typeof(snp[:gene]))
            #println(typeof(snp[:rsid]))
            if !in(snp[:gene], keys(gene_snps))
                gene_snps[snp[:gene]] = [snp[:rsid]]
            else
                append!(gene_snps[snp[:gene]],[snp[:rsid]])
            end
        end
        gene_info = DataFrame(gene_id=String[],num_ind_any=Int64[],max_ind_per_snp=Int64[],num_snps_any=Int64[],max_snp_per_ind=Int64[])
        for gene in keys(gene_snps)
            ind_set = Set{String}()
            max = 0
            missing_snps = length(gene_snps[gene])
            for snp in gene_snps[gene]
                if !in(snp,keys(snp_dict)) # number of SNPs missing in at least one person for each model
                    missing_snps -= 1
                    continue
                end
                for ind in snp_dict[snp] # number of people missing a SNP in each model
                    push!(ind_set,ind)
                end
                if length(snp_dict[snp]) > max # max number of people missing a single SNP per model
                    max = length(snp_dict[snp])
                end        
        # max number of missing SNPs in a single person for each model
            end
        end
    end
end

function main()
    parsed_args = parse_commandline()
    if parsed_args["performance"]
        performance(parsed_args["databases"])
    end
    if parsed_args["missingness"]
        missingSNPs(parsed_args["databases"],parsed_args["genomes"])
    end
end

main()
