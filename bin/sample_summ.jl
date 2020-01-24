# quick script to get sample summary by num_snps threshold.
# adapted from functions in snp_choice_analyses.jl

using CSV
using DataFrames
using Plots
using Seaborn
using StatsBase

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

function mapSampleIDs(df::DataFrames.DataFrame,map_file::String,info_file::String)
    map_ids = CSV.read(map_file)
    df = join(df,map_ids,on=:master_id => :reich_id,kind=:inner)
    df = join(df,CSV.read(info_file),kind=:inner,on=:colbran_id => :sample,makeunique=true)
    return df
end

function sampleSummary()
    anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    map_file =  "data/reich_ancient_humans/mapped_sample_ids.csv"
    info_file = "data/aDNA_sources.csv"
    anno = parseAnnoFile(anno_file)
    anno = anno[anno[:assessment] .!= "no",:]
    println("Number of Samples at Start: $(nrow(anno))")
    anno = mapSampleIDs(anno,map_file,info_file)
    println("Number of Samples successfully mapped: $(nrow(anno))")
    anno[anno[:lifestyle] .== "nomadic/pastoral",:lifestyle] = "Past."
    anno[anno[:lifestyle] .== "pastoral",:lifestyle] = "Past."
    anno[anno[:lifestyle] .== "Agriculture",:lifestyle] = "Agri."
    anno[anno[:lifestyle] .== "Hunter-gatherer",:lifestyle] = "HG"

    # for ones with >777677 SNPs called (3rd quartile)
    anno = anno[anno[:num_snps].>777677,:]
    println("with >777677:")
    describe(anno[:num_snps])
    describe(anno[:date])
    f, axes = Seaborn.subplots(2, 3, figsize=(10,7), sharey="all")
    x = y = 2
    plot_num = 1
    println("\nby continent:")
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
    Seaborn.savefig("lifestyle_age_by_continent_3rdq.pdf")
    clf()
    CSV.write("reich_individuals_3rdq.txt",anno[[:ind_id,:master_id,:paper,:date,:group_id,:country,:sex,:coverage,:num_snps,:continent,:colbran_id,:type,:location,:population,:lifestyle,:latitude,:longitude]];delim='\t')

end

sampleSummary()
