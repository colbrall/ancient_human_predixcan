# gene_or.jl
#
# 04/13/2020
# @author Laura Colbran
# calculates OR of gene categories using output from compare_by_group.jl
#

using DataFrames
using CSV
using ArgParse
using GZip
using HypothesisTests
using StatsPlots

#set plotting defaults
default(color=:black,leg=false,grid=false,fontfamily="arial")

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--path","-p"
            help = "path to file with pos/neg calls"
            arg_type = String
            required=true
        "--column","-c"
            help = "column of file that has true/false calls. assumes gene/tissue are columns 1/2."
            arg_type = Int64
            required=true
        "--targets","-t"
            help = "file of genes in the target set with column called 'ensembl_id'"
            arg_type = String
        "--tissue"
            help = "if you want to do odds ratio by tissue, rather than overall"
            action=:store_true
        "--summarize"
            help = "if you want to plot summary distributions"
            action=:store_true
    end
    return parse_args(s)
end

function overallOR(file_path::String,col::Int64,targ_file::String)
    # calulate Odds Ratio of Target for each set vs. everything else
    genes = CSV.read(file_path;delim='\t',allowmissing=:none)[:,[1,col]]
    genes = by(genes,:gene,names(genes)[2]=> x-> any(x))
    targets = CSV.read(targ_file;delim='\t',allowmissing=:none)[:,[:ensembl_id]]
    genes = join(genes, targets, on=:gene=>:ensembl_id,indicator = :source,kind=:left)
    sig_targ = nrow(genes[(genes[:source].=="both").&(genes[names(genes)[2]].==true),:])
    sig_nontarg = nrow(genes[(genes[:source].=="left_only").&(genes[names(genes)[2]].==true),:])
    nonsig_targ = nrow(genes[(genes[:source].=="both").&(genes[names(genes)[2]].==false),:])
    nonsig_nontarg = nrow(genes[(genes[:source].=="left_only").&(genes[names(genes)[2]].==false),:])
    println("\ncontingency table is structured:")
    println("    Sig. Target | Sig. not-Target")
    println("Not Sig. Target | Not Sig. not-Target\n")
    println(FisherExactTest(sig_targ,sig_nontarg,nonsig_targ,nonsig_nontarg))
end

function tissueOR(file_path::String,col::Int64,targ_file::String)
    # calulate Odds Ratio of Target for each set vs. everything else
    genes = CSV.read(file_path;delim='\t',allowmissing=:none)[:,[1,2,col]]
    targets = CSV.read(targ_file;delim='\t',allowmissing=:none)[:,[:ensembl_id]]
    genes = join(genes, targets, on=:gene=>:ensembl_id,indicator = :source,kind=:left)
    for tiss in unique(genes[:tiss])
        temp = genes[genes[:tiss] .== tiss,:]
        sig_targ = nrow(temp[(temp[:source].=="both").&(temp[names(temp)[3]].==true),:])
        sig_nontarg = nrow(temp[(temp[:source].=="left_only").&(temp[names(temp)[3]].==true),:])
        nonsig_targ = nrow(temp[(temp[:source].=="both").&(temp[names(temp)[3]].==false),:])
        nonsig_nontarg = nrow(temp[(temp[:source].=="left_only").&(temp[names(temp)[3]].==false),:])
        println(tiss)
        println("\ncontingency table is structured:")
        println("    Sig. Target | Sig. not-Target")
        println("Not Sig. Target | Not Sig. not-Target\n")
        println(FisherExactTest(sig_targ,sig_nontarg,nonsig_targ,nonsig_nontarg))
    end
end

function summarize(file_path::String,col::Int64)
    all = CSV.read(file_path;delim='\t',allowmissing=:none)[:,[1,2,col]]
    all[:count] = 1
    plot_df = by(all[all[3].==true,:],:gene,:count=>sum)
    describe(plot_df[:count_sum])
    gene_plot = histogram(plot_df[:count_sum],xlabel="Number of Tissues",ylabel="Number of Genes",bins=100,margin=10Plots.mm)
    Plots.savefig(gene_plot,"gene_numtiss_dist.pdf")
    plot_df = by(all[all[3].==true,:],:tiss,:count=>sum)
    plot_df[:tot_sum] = 0
    for i in 1:nrow(plot_df)
        plot_df[i,:tot_sum] = nrow(all[all[2].==plot_df[i,1],:])
    end
    plot_df[:prop] = plot_df[:count_sum] ./ plot_df[:tot_sum]
    sort!(plot_df,names(plot_df)[1])
    println(plot_df[1])
    describe(plot_df[:prop])
    tiss_plot = bar(plot_df[1],plot_df[:prop],xlabel="Tissues",ylabel="Prop. Sig. Genes",margin=10Plots.mm)
    Plots.savefig(tiss_plot,"tiss_numgene.pdf")
end

function main()
    parsed_args = parseCommandLine()
    if parsed_args["tissue"]
        tissueOR(parsed_args["path"],parsed_args["column"],parsed_args["targets"])
    elseif parsed_args["summarize"]
        summarize(parsed_args["path"],parsed_args["column"])
    else
        overallOR(parsed_args["path"],parsed_args["column"],parsed_args["targets"])
    end
end

main()
