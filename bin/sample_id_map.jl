# sample_map.jl
# @author Laura Colbran, May 2019
#
# makes the sample map from my aDNA_sources.CSV to reich_lab IDs in their compilation
#
# Julia 1.1

using CSV
using DataFrames

EXCLUDE = ["Ancestor.REF","Href.REF","Chimp.REF", "Gorilla.REF","Altai_snpAD.DG",
            "Denisova_snpAD.DG","Altai_published.DG","Denisova_published.DG",
            "Denisova11.SG","Vindija_snpAD.DG","VindijaG1_final.SG","Goyet_final_provisional.SG",
            "Les_Cottes_final_provisional.SG","Mezmaiskaya1_final_provisional.SG","VindijaG1_final_provisional.SG"]

# parses Reich Lab sample annotation file
function parseAnnoFile(anno_file::String)
    inds = CSV.read(anno_file; delim='\t',allowmissing=:none)
    inds = inds[[2,3,4,7,9,16],]
#    inds = inds[[1,2,3,7,8,9],]
    names!(inds,[:ind_id,:master_id,:skeleton_id,:paper,:date,:type])
#   remove present-day samples
    println(nrow(inds))
#    inds = inds[inds[:date] .!= "0",:]
    inds = inds[inds[:date] .!= 0,:]
    println(nrow(inds))

#   remove archaic hominins, reference genomes
    inds[:keep] = ".."
    for i in 1:nrow(inds)
        if inds[i,:ind_id] in EXCLUDE
            inds[i,:keep] = "NO"
        end
    end
    inds = inds[inds[:keep] .!== "NO",:]
    deletecols!(inds, :keep)
    return inds
end

function main()
    anno_file = "../../data/ancient_dna/reich_compilation/v42/v42.4.1240K.anno"
    source_file = "data/aDNA_sources.csv"
    anno = parseAnnoFile(anno_file)
    source = CSV.read(source_file)
#   remove leading/trailing whitespace
    for i in 1:nrow(source)
        source[i,:sample] = strip(source[i,:sample])
    end
#   all the simple joins; my sample id == master_id
    id_map = DataFrames.DataFrame(colbran_id = join(source,anno,on=:sample => :master_id,kind=:inner,makeunique=true)[:sample])
    id_map[:reich_id] = id_map[:colbran_id]
#   grab the IDs the simple way  didn't work for
    problems = join(anno,source,on=:master_id => :sample,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   joins where my id matched a skeleton id, split by a comma
    problems[:sk_short] = ""
    for i in 1:nrow(problems)
        problems[i,:sk_short] = split(problems[i,:skeleton_id],", ")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   joins where my id matched a skeleton id, split by a comma, different order than previous
    for i in 1:nrow(problems)
        problems[i,:sk_short] = split(problems[i,:skeleton_id],", ")[end]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:sample],
                  reich_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   ind_id was a modified version of skeleton id: "_" different position
    for i in 1:nrow(problems)
        problems[i,:sk_short] = split(problems[i,:skeleton_id],"_")[end]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:sample],
                  reich_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   sample id == :ind_id
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :ind_id,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(source,problems,on=:sample => :ind_id,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   ind_id was a modified version of sample id: "."
    problems[:ind_short] = ""
    for i in 1:nrow(problems)
        problems[i,:ind_short] = split(problems[i,:ind_id],".")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :ind_short,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(source,problems,on=:sample => :ind_short,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   ind_id was a modified version of sample id: "_"
    for i in 1:nrow(problems)
        problems[i,:ind_short] = split(problems[i,:ind_id],"_")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :ind_short,kind=:inner,makeunique=true)[:sample],
                  reich_id = join(source,problems,on=:sample => :ind_short,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   sample id was a modified version of :master_id - "*"
    source[:id_short] = ""
    for i in 1:nrow(source)
        source[i,:id_short] = split(source[i,:sample],"*")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:id_short => :master_id,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(source,problems,on=:id_short => :master_id,kind=:inner,makeunique=true)[:id_short]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")

#   sample id was a modified version of :master_id -  "-"
    source[:id_short] = ""
    for i in 1:nrow(source)
        source[i,:id_short] = split(source[i,:sample],"-")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:id_short => :master_id,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(source,problems,on=:id_short => :master_id,kind=:inner,makeunique=true)[:id_short]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")

#   hand-matched
    source[:spot_check] = ""
    for i in 1:nrow(source)
        if source[i,:sample] == "SZ27B"
            source[i,:spot_check] = "SZ27"
        elseif source[i,:sample] == "RISE508**"
            source[i,:spot_check] = "RISE507.508.merge"
        elseif source[i,:sample] == "GB1_Eneo"
            source[i,:spot_check] = "GB"
        elseif source[i,:sample] == "Chan_Meso"
            source[i,:spot_check] = "Chan"
        elseif source[i,:sample] == "Canes_Meso"
            source[i,:spot_check] = "Canes"
        elseif source[i,:sample] == "JomonB"
            source[i,:spot_check] = "Sanganji131464.SG"
        elseif source[i,:sample] == "JomonA"
            source[i,:spot_check] = "Sanganji131421-3"
        elseif source[i,:sample] == "Bot15"
            source[i,:spot_check] = "Botocudo15.SG"
        elseif source[i,:sample] == "3DRIF-16"
            source[i,:spot_check] = "3DT16"
        elseif source[i,:sample] == "3DRIF-26"
            source[i,:spot_check] = "3DT26"
        elseif source[i,:sample] == "6DRIF-18"
            source[i,:spot_check] = "6DT18"
        elseif source[i,:sample] == "6DRIF-21"
            source[i,:spot_check] = "6DT21"
        elseif source[i,:sample] == "6DRIF-22"
            source[i,:spot_check] = "6DT22"
        elseif source[i,:sample] == "6DRIF-23"
            source[i,:spot_check] = "6DT23"
        elseif source[i,:sample] == "6DRIF-3"
            source[i,:spot_check] = "6DT3"
        elseif source[i,:sample] == "NO3432"
            source[i,:spot_check] = "I17268"
        elseif source[i,:sample] == "I2656_d"
            source[i,:spot_check] = "I2656"
        elseif source[i,:sample] == "I8140_d"
            source[i,:spot_check] = "I8140"
        elseif source[i,:sample] == "I9058_d"
            source[i,:spot_check] = "I9058"
        elseif source[i,:sample] == "939"
            source[i,:spot_check] = "XVII-B-939"
        elseif source[i,:sample] == "Cinchorro"
            source[i,:spot_check] = "Chinchorroi15" #N.B. the sex differs between mine and reich. date/cov/loc agree, and it's the only one in both sets, so one of us must be wrong.
        elseif source[i,:sample] == "Saqqaq"
            source[i,:spot_check] = "Inuk"
        elseif source[i,:sample] == "Kostenki14_wgs"
            source[i,:spot_check] = "Kostenki14"
        elseif source[i,:sample] == "I2966_1240k" #for sake of lifestyle info
            source[i,:spot_check] = "I2966_all"
        elseif source[i,:sample] == "I2967_1240k" #for sake of lifestyle info
            source[i,:spot_check] = "I2967_all"
        elseif source[i,:sample] == "Be9_1240k" #for sake of lifestyle info
            source[i,:spot_check] = "I0562"
        elseif source[i,:sample] == "Ze6"
            source[i,:spot_check] = "Ze6b"
        elseif source[i,:sample] == "A10"
            source[i,:spot_check] = "I0577"
        elseif source[i,:sample] == "Be11"
            source[i,:spot_check] = "I0563"
        elseif source[i,:sample] == "A17"
            source[i,:spot_check] = "I0576"
        elseif source[i,:sample] == "PR9"
            source[i,:spot_check] = "I0574"
        elseif source[i,:sample] == "PR3"
            source[i,:spot_check] = "I0575"
        elseif source[i,:sample] == "Alh"
            source[i,:spot_check] = "Alh3a"
        elseif source[i,:sample] == "BAL001" #for sake of lifestyle info
            source[i,:spot_check] = "BAL051"
        elseif source[i,:sample] == "SIJ002.A0101"
            source[i,:spot_check] = "I11131"
        elseif source[i,:sample] == "SIJ003.A0101"
            source[i,:spot_check] = "I11132"
        elseif source[i,:sample] == "SIJ001.A01"
            source[i,:spot_check] = "I11133"
        elseif source[i,:sample] == "VEK007.A0101"
            source[i,:spot_check] = "VEK007"
        elseif source[i,:sample] == "JomonA"
            source[i,:spot_check] = "Sanganji_131421-3_A1"
        elseif source[i,:sample] == "JomonB"
            source[i,:spot_check] = "Sanganji_131464_B"
        end
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:spot_check => :master_id,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(source,problems,on=:spot_check => :master_id,kind=:inner,makeunique=true)[:spot_check]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")

    println(sort(unique(problems[:paper])))
    println(problems[problems[:paper].=="SchroederPNAS2019",:])


    CSV.write("data/reich_ancient_humans/mapped_sample_ids.csv", unique(id_map);delim=',')
end

main()
