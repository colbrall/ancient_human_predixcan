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

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--databases","-d"
            nargs=2
            help = "paths to the 1 or 2 model databases"
            arg_type = String
            required = true
        "--performance"
            help = "to run comparison between certain models' performances and SNP counts across different training instances"
            action = :store_true
        "--snps"
            help = "to compare specific SNPs and their weights"
            action = :store_true
        "--genomes","-g"
            help = "path to directory containing genomes"
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
	println("Original models first, then mine:")
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

function main()
    parsed_args = parse_commandline()
    if parsed_args["performance"]
        performance(parsed_args["databases"])
    end
end

main()
