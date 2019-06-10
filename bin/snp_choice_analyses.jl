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
using Plots
using Seaborn
using StatsBase
using HypothesisTests

#set plotting defaults
default(color=:black,leg=false,grid=false,fontfamily="arial")

Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))



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
        "--summary"
            help = "to calculate summary stats for individuals in annotation file"
            action = :store_true

        "--anno_file","-a"
            help = "path to sample annotation file"
            arg_type = String

        "--snp_coverage"
            help = "for each SNP, calculates coverage metric across samples, then summmarises and outputs"
            action = :store_true

        "--geno_dir","-g"
            help = "directory containing dosage files for samples. required for snp_coverage."
            arg_type = String

        "--ind_file","-i"
            help = "fam file corresponding to the individuals in your dosage file. order matters."
            arg_type = String

        "--snp_summary"
            help = "if you want to run summary stats on SNP coverage"
            action = :store_true

        "--snp_file","-s"
            help = "bed file of SNPs containing their coverage rate and genotype count. Output of snpCoverage() or subset of it"
            arg_type = String

        "--sample_summary"
            help = "if you want to run summary stats on samples based on a given set of SNPs of interest and a threshold for missing data"
            action = :store_true

        "--map_file","-m"
            help = "csv file containing ID mapping. output of sample_id_map.jl."
            arg_type = String

        "--info_file","-n"
            help = "file containing further sample info. aDNA_sources.csv usually."
            arg_type = String

        "--threshold","-t"
            help = "highest proportion of missing SNPs allowed for a sample to be considered if calling --sample_summary. Lowest allowable coverage rate if --pick_snps."
            arg_type = Float64
            default = 1.0

        "--pick_snps"
            help = "to extract SNPs passing a certain coverage rate"
            action = :store_true
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

# calculates coverage rate for a SNP
function coverageRate(genotypes::Array{SubString{String},1},num_snps::Array{Int64,1})
    rate = sum(num_snps[findall(!isequal("NA"),genotypes)])
    return rate
end

# counts number of samples with calls for a SNP
function countSamples(genotypes::Array{SubString{String}, 1})
    ct = length(findall(!isequal("NA"),genotypes))
    return ct
end

# reads in .fam file and returns ids
function indArray(ind_file::String)
    ids = Array{String,1}()
    for line in readlines(open(ind_file,"r"))
        push!(ids,String(split(line,' ')[2]))
    end
    return ids
end

# boolean whether the genotyped snps in a sample is a high enough proportion to pass threshold
function passThresh(count::Int64,total::Int64,thresh::Float64)
    if count/total < 1-thresh
        return false #if proportion called is fewer than necessary
    else
        return true
    end
end

# returns array of :ind_ids passing a threshold of proportion of SNPs covered
function filterThresh(snps::String,geno_dir::String,ind_file::String,thresh::Float64)
    snp_ids = CSV.read(snps;delim='\t')[:rsid]
    inds = indArray(ind_file)
    counts = zeros(Int64,length(inds),1)
    for item in filter(x -> occursin(r"(?i)\.txt\.gz", x), readdir(geno_dir))
        println(item)
        for line in readlines(GZip.open("$geno_dir$item"))
            if !(String(split(line," ")[2]) in snp_ids) continue end #skip snps not of interest
            genotypes = split(line,' ')[7:end]
            counts[findall(!isequal("NA"),genotypes)] .+= 1
        end
    end
    return inds[findall(i -> passThresh(i,length(snp_ids),thresh),counts)]
end

# maps sample ids in my Spreadsheet O' Knowledge to those in annotation file
# also joins SoK to annotation file by that map
function mapSampleIDs(df::DataFrames.DataFrame,map_file::String,info_file::String)
    map_ids = CSV.read(map_file)
    df = join(df,map_ids,on=:master_id => :reich_id,kind=:inner)
    df = join(df,CSV.read(info_file),kind=:inner,on=:colbran_id => :sample,makeunique=true)
    return df
end

#summary stats of samples
function summaryInds(anno_file::String)
    anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    anno = parseAnnoFile(anno_file)
    println("Total Ancient Human Samples: $(nrow(anno))")
    println("Number Reich Lab considered 'PASS': $(nrow(anno[anno[:assessment].== "yes",:]))")
    println("Number they considered seriously questionable: $(nrow(anno[anno[:assessment] .== "no",:]))")
    anno = anno[anno[:assessment] .!= "no",:]
    println("Number of unique individuals: $(length(unique(anno[:master_id])))")
    println("\nNumber with date information: $(nrow(anno[anno[:date] .> 0,:]))")
    println("Number with coverage information: $(nrow(anno[anno[:coverage] .> 0,:]))\n")
    describe(anno[:date])
    age_plot = histogram(anno[:date], nbins=100,xlabel="Age (years BP)",ylabel="Number of Samples")

    describe(anno[:coverage])
    cov_plot = histogram(anno[:coverage],nbins=100,xlabel="Coverage",ylabel="Number of Samples")

    describe(anno[:num_snps])
    snp_plot = histogram(anno[:num_snps], nbins=100, ylabel="Number of Samples",xlabel= "Number Autosomal SNPs Called")

    println("\nRho: $(corspearman(anno[:coverage],anno[:num_snps]))")
    println(OneSampleZTest(atanh(corspearman(anno[:coverage],anno[:num_snps])),1, nrow(anno)))
    println("Number with >30k SNPs: $(nrow(anno[anno[:num_snps] .> 30000,:]))\n")

    for cont in unique(anno[:continent])
        println("Num. samples in $cont: $(nrow(anno[anno[:continent] .== cont,:]))")
    end
    for type in unique(anno[:type])
        println("Num. $type samples: $(nrow(anno[anno[:type] .== type,:]))")
    end

    CSV.write("data/reich_ancient_humans/ancient_human_ids.txt",DataFrames.DataFrame(id = anno[:ind_id]); delim='\t')
    Plots.savefig(cov_plot, "sample_coverage_distribution.pdf")
    Plots.savefig(age_plot,"sample_date_distribution.pdf")
    Plots.savefig(snp_plot, "sample_numsnps_distribution.pdf")
end

# calculates and writes out coverage rate: Sum(Num_Snps) across individuals that
# have a call for the site where Num_Snps is the overall coverage of each individual
function snpCoverage(geno_dir::String,anno_file::String,ind_file::String)
    #anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    #ind_file = "data/reich_ancient_humans/kept_ancient_humans.fam"
    snps = DataFrames.DataFrame(chr = String[],pos = Int64[],rsid = String[],
                                        coverage_rate = Int64[],gt_count = Int64[])
    samples = parseAnnoFile(anno_file)
    samples = join(deletecols!(CSV.read(ind_file;delim=' ',
                    header=[:b,:ind_id,:x,:y,:a,:z],allowmissing=:none),[:a,:b,:x,:y,:z]),
                    samples,kind=:inner,on=:ind_id) #filter to just inds of interest
    samples[:num_snps]
    for item in filter(x -> occursin(r"(?i)\.txt\.gz", x), readdir(geno_dir))
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
    dist_plot = histogram(snps[:coverage_rate], nbins = 100,xlabel = "Coverage Rate",ylabel = "Number of SNPs")
    Plots.savefig(dist_plot, "1240k_snps_coverage_rate.pdf")
    println("\nGenotype Count:")
    describe(snps[:gt_count])
    dist_plot = histogram(snps[:gt_count], nbins = 100,xlabel = "Number of Samples with Genotype",ylabel = "Number of SNPs")
    Plots.savefig(dist_plot, "1240k_snps_genotype_count.pdf")

    # dist_plot = histogram(snps[snp[:type].== "1240K",:coverage_rate], nbins = 100,xlabel = "Coverage Rate",ylabel = "Number of SNPs",color=:red,alpha=0.3)
    # dist_plot = histogram!(snps[snp[:type].!= "1240K",:coverage_rate], nbins = 100,xlabel = "Coverage Rate",ylabel = "Number of SNPs",color=:blue,alpha=0.3)
    # Plots.savefig(dist_plot, "1240k_snps_coverage_rate.pdf")
end

function sampleSummary(anno_file::String,snp_file::String, geno_dir::String,ind_file::String,threshold::Float64,info_file::String,map_file::String)
    anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    map_file =  "data/reich_ancient_humans/mapped_sample_ids.csv"
    info_file = "data/aDNA_sources.csv"
    anno = parseAnnoFile(anno_file)
    anno = anno[anno[:assessment] .!= "no",:]
    println("Number of Samples at Start: $(nrow(anno))")
    passed = filterThresh(snp_file,geno_dir,ind_file,threshold)
    anno = anno[[in(id,passed) for id in anno[:ind_id]],:]
    println("Number of Samples passing coverage threshold: $(nrow(anno))")
    anno = mapSampleIDs(anno,map_file,info_file)
    println("Number of Samples successfully mapped: $(nrow(anno))")
    anno[anno[:lifestyle] .== "nomadic/pastoral",:lifestyle] = "Past."
    anno[anno[:lifestyle] .== "pastoral",:lifestyle] = "Past."
    anno[anno[:lifestyle] .== "Agriculture",:lifestyle] = "Agri."
    anno[anno[:lifestyle] .== "Hunter-gatherer",:lifestyle] = "HG"
    f, axes = Seaborn.subplots(2, 3, figsize=(10,7), sharey="all")
    x = y = 2
    plot_num = 1
    for loc in unique(anno[:continent])
        println("\n$loc")
        for life in unique(anno[:lifestyle])
            println("Num $life: $(count(i->(i==life),anno[anno[:continent].==loc,:lifestyle]))")
        end
        describe(anno[anno[:continent].==loc,:date])
        life_plot = swarmplot(anno[anno[:continent] .== loc,:lifestyle],anno[anno[:continent] .== loc,:date],
                            order=["HG","Past.","Agri.","NA"],ax=axes[x, y])
        life_plot.set_title("$(loc) (N = $(nrow(anno[anno[:continent] .== loc,:])))")
        life_plot.set_ylabel("Age (Years BP)")
        plot_num += 1
        x = plot_num%2 + 1
        y = plot_num%3 + 1
    end
    Seaborn.savefig("lifestyle_age_by_continent.pdf")
    clf()
end

# extracts and outputs SNPs that pass a certain coverage rate threshold.
function extractSNPs(snp_file::String,thresh::Float64)
    snps = CSV.read(snp_file;delim='\t')
    CSV.write("$(splitext(snp_file)[1])_filtered_$(thresh).bed",snps[snps[:coverage_rate] .> thresh,:];delim='\t')
end

function main()
    args = parse_commandline()
    if args["summary"]
        summaryInds(args["anno_file"])
    end
    if args["snp_coverage"]
        snpCoverage(args["geno_dir"],args["anno_file"],args["ind_file"])
    end
    if args["snp_summary"]
        snpSummary(args["snp_file"])
    end
    if args["sample_summary"]
        sampleSummary(args["anno_file"],args["snp_file"],args["geno_dir"],args["ind_file"],
                        args["threshold"],args["info_file"],args["map_file"])
    end
    if args["pick_snps"]
        extractSNPs(args["snp_file"],args["threshold"])
    end
end

main()
