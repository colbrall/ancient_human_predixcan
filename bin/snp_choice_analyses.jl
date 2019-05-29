# snp_choice_analyses.jl
#
# @author Laura Colbran, 5/22/2019
#
# tests run to decide which set of SNPs to retrain PrediXcan models on
# also more general sample exploration
# julia 1.1

using ArgParse
using CSV
using DataFrames
using GZip
using StatsPlots
using StatsBase
using HypothesisTests
pyplot(grid = false,leg = false)

EUROPE = ["Austria","Belgium","Bulgaria","Croatia","Czech Republic","Denmark",
            "Estonia","Finland","France","","Germany","Great Britain",
            "Greece","Hungary","Iceland","Ireland","Italy","Kosovo","Latvia",
            "Liechtenstein","Lithuania","Luxembourg","Macedonia","Malta","Moldova",
            "Monaco","Montenegro","Netherlands","Norway","Poland","Portugal",
            "Romania","San Marino","Serbia","Slovakia","Slovenia","Spain","Sweden",
            "Switzerland","Ukraine","Vatican"] #countries in europe in annotation file
ASIA = ["Armenia", "Cambodia", "China","Georgia","India","Iran","Israel","Japan",
        "Jordan", "Kazakhstan", "Kyrgyzstan", "Laos", "Lebanon","Mongolia",
        "Myanmar", "Nepal", "Russia","Thailand", "Turkey", "Turkmenistan","Vietnam"]
        #countries in asia in annotation file
N_AMERICA = ["Bahamas", "Belize", "Canada","Greenland","Mexico", "USA"] #countries in north america in annotation file
S_AMERICA = ["Argentina","Brazil","Chile", "Peru"] #countries in south america in annotation file
OCEANIA = ["French Polynesia", "Indonesia","Malaysia", "Philippines", "Solomon Islands",
            "Vanuatu"] #countries in oceania in annotation file
AFRICA = ["Canary Islands", "Egypt", "Ethiopia","Kenya","Malawi", "Morocco",
        "South Africa", "Tanzania","Tonga"] #countries in africa in annotation file
EXCLUDE = ["Ancestor.REF","Href.REF","Chimp.REF", "Gorilla.REF","Altai_snpAD.DG",
            "Denisova_snpAD.DG","Altai_published.DG","Denisova_published.DG",
            "Denisova11.SG","Vindija_snpAD.DG","VindijaG1_final.SG"] #non-human ids in annotation file

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--anno_file","-a"
            help = "path to sample annotation file"
            arg_type = String

        "--summary"
            help = "to calculate summary stats for individuals in annotation file"
            action = :store_true

        "--snp_coverage"
            help = "for each SNP, calculates coverage metric across samples, then summmarises and outputs"
            action = :store_true

        "--geno_file","-g"
            help = "dosage file containing genotypes for samples. required for snp_coverage."
            arg_type = String

        "--ind_file","-i"
            help = "fam file corresponding to the individuals in your dosage file. order matters."
            arg_type = String

        "--snp_summary"
            help = "if you want to run summary stats on SNP coverage"
            action = :store_true

        "--snp_file","-s"
            help = "bed file of SNPs containing their coverage rate and genotype count. Output of snpCoverage()"
            arg_type = String
    end
    return parse_args(s)
end

# parses Reich Lab sample annotation file
function parseAnnoFile(anno_file::String)
    inds = CSV.read(anno_file; delim='\t',allowmissing=:none)
    inds = inds[[1,2,7,8,9,11,13,16,20,21,24],]
    names!(inds,[:ind_id,:master_id,:type,:paper,:date,:group_id,:country,:sex,:coverage,:num_snps,:assessment])
    inds = inds[inds[:date] .!= "0",:] #remove present-day samples
    inds[:continent] = ".."
    inds[:ind_date] = 0
    inds[:flt_cov] = 0.0
    for i in 1:nrow(inds)
        if inds[i,:date] == ".."
            inds[i,:ind_date] = -1
        else
            inds[i,:ind_date] = parse(Int,inds[i,:date])
        end

        if inds[i,:coverage] == ".."
            inds[i,:flt_cov] = -1.0
        else
            inds[i,:flt_cov] = parse(Float64,inds[i,:coverage])
        end

        if inds[i,:country] in EUROPE
            inds[i,:continent] = "Europe"
        elseif inds[i,:country] in ASIA
            inds[i,:continent] = "Asia"
        elseif inds[i,:country] in N_AMERICA
            inds[i,:continent] = "N_America"
        elseif inds[i,:country] in S_AMERICA
            inds[i,:continent] = "S_America"
        elseif inds[i,:country] in AFRICA
            inds[i,:continent] = "Africa"
        elseif inds[i,:country] in OCEANIA
            inds[i,:continent] = "Oceania"
        end

        if inds[i,:ind_id] in EXCLUDE
            inds[i,:continent] = "NON-HUMAN"
        end

        if startswith(inds[i,:assessment],"PASS")
            inds[i,:assessment] = "yes"
        elseif startswith(inds[i,:assessment],"QUESTIONABLE_CRITICAL")
            inds[i,:assessment] = "no"
        elseif startswith(inds[i,:assessment],"QUESTIONABLE")
            inds[i,:assessment] = "maybe"
        end
    end
    inds[:date] = inds[:ind_date]
    inds[:coverage] = inds[:flt_cov]
    deletecols!(inds,[:ind_date,:flt_cov])
    inds = inds[inds[:continent] .!== "NON-HUMAN",:] #remove archaic hominins, reference genomes
    return inds
end

#summary stats of samples
function summaryInds(anno_file::String)
    #anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    anno = parseAnnoFile(anno_file)
    println("Total Ancient Human Samples: $(nrow(anno))")
    println("Number Reich Lab considered 'PASS': $(nrow(anno[anno[:assessment].== "yes",:]))")
    println("Number they considered seriously questionable: $(nrow(anno[anno[:assessment] .== "no",:]))")
    anno = anno[anno[:assessment] .!= "no",:]
    println("Number of unique individuals: $(length(unique(anno[:master_id])))")
    CSV.write("data/reich_ancient_humans/ancient_human_ids.txt",DataFrames.DataFrame(id = anno[:ind_id]); delim='\t')
    println("\nNumber with date information: $(nrow(anno[anno[:date] .> 0,:]))")
    println("Number with coverage information: $(nrow(anno[anno[:coverage] .> 0,:]))\n")
    describe(anno[:date])
    dist_plot = histogram(anno[:date], nbins = 100,
                    xlabel = "Date (years BP)",
                    ylabel = "Count",
                    color = [:black])
    savefig(dist_plot, "sample_date_distribution.pdf")
    describe(anno[:coverage])
    dist_plot = histogram(anno[:coverage], nbins = 100,
                    xlabel = "Coverage",
                    ylabel = "Count",
                    color = [:black])
    savefig(dist_plot, "sample_coverage_distribution.pdf")
    describe(anno[:num_snps])
    dist_plot = histogram(anno[:num_snps], nbins = 100,
                    xlabel = "Number Autosomal SNPs",
                    ylabel = "Count",
                    color = [:black])
    savefig(dist_plot, "sample_numsnps_distribution.pdf")

    println("\nRho: $(corspearman(anno[:coverage],anno[:num_snps]))")
    println(OneSampleZTest(atanh(corspearman(anno[:coverage],anno[:num_snps])),1, nrow(anno)))
    println("Number with >30k SNPs: $(nrow(anno[anno[:num_snps] .> 30000,:]))\n")

    for cont in unique(anno[:continent])
        println("Num. samples in $cont: $(nrow(anno[anno[:continent] .== cont,:]))")
    end
end

# calculates coverage rate for a SNP
function coverageRate(genotypes::Array{SubString{String},1},num_snps::Array{Int64,1})
    #println(genotypes)
    rate = sum(num_snps[findall(!isequal("NA"),genotypes)])
    return rate
end

# counts number of samples with calls for a SNPs
function countSamples(genotypes::Array{SubString{String}, 1})
    ct = length(findall(!isequal("NA"),genotypes))
    return ct
end

# calculates and writes out coverage rate: Sum(Num_Snps) across individuals that
# have a call for the site where Num_Snps is the overall coverage of each individual
function snpCoverage(geno_dir::String,anno_file::String,ind_file::String)
    #anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    #ind_file = "data/reich_ancient_humans/kept_ancient_humans.fam"
    snps = DataFrames.DataFrame(chr = String[],pos = Int64[],rsid = String[],coverage_rate = Int64[],gt_count = Int64[])
    samples = parseAnnoFile(anno_file)
    samples = join(deletecols!(CSV.read(ind_file;delim=' ',header=[:b,:ind_id,:x,:y,:a,:z],allowmissing=:none),[:a,:b,:x,:y,:z]),
                    samples,kind=:inner,on=:ind_id) #filter to just inds of interest
    samples[:num_snps]
    for item in readdir(geno_dir)
        if !endswith(item,"txt.gz") continue end
        println(item)
        for line in readlines(GZip.open("$geno_dir$item"))
            push!(snps,["chr$(split(line,' ')[1])",parse(Int,split(line,' ')[3]),split(line,' ')[2],
                        coverageRate(split(line,' ')[7:end],samples[:num_snps]),
                        countSamples(split(line,' ')[7:end])])
        end
    end
    snps[:start] = snps[:pos] .- 1
    CSV.write("snp_coverage_rate.bed",snps[[1,6,2,3,4,5]];delim='\t')
end

# summary stats and plot for snp coverage (output from snpCoverage())
function snpSummary(snp_file::String)
    snps = CSV.read(snp_file;delim='\t')
    println("Coverage Rate:")
    describe(snps[:coverage_rate])
    dist_plot = histogram(snps[:coverage_rate], nbins = 100,
                    xlabel = "Coverage Rate",
                    ylabel = "Number of SNPs",
                    color = [:black])
    savefig(dist_plot, "1240k_snps_coverage_rate.pdf")
    println("\nGenotype Count:")
    describe(snps[:gt_count])
    dist_plot = histogram(snps[:gt_count], nbins = 100,
                    xlabel = "Number of Samples with Genotype",
                    ylabel = "Number of SNPs",
                    color = [:black])
    savefig(dist_plot, "1240k_snps_genotype_count.pdf")
end

# maps sample ids in my Spreadsheet O' Knowledge to those in annotation file
function mapSampleIDs()
    println("Woops not implemented yet!")
    # anno id of interest is master_id
end

function sampleSummary()
    println("Woops not implemented yet!")
    # read in list of samples with enough coverage in SNPs we care about
    # join to anno_file for info
    # map ids to ids in aDNA_sources
    # join to aDNA_sources
    # make plots corinne described-- split by location, with lifestyle and yearBP
end

function main()
    args = parse_commandline()
    if args["summary"]
        summaryInds(args["anno_file"])
    end
    if args["snp_coverage"]
        snpCoverage(args["geno_file"],args["anno_file"],args["ind_file"])
    end
    if args["snp_summary"]
        snpSummary(args["snp_file"])
    end
end

main()
