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
            "Denisova11.SG","Vindija_snpAD.DG","VindijaG1_final.SG"]

# parses Reich Lab sample annotation file
function parseAnnoFile(anno_file::String)
    inds = CSV.read(anno_file; delim='\t',allowmissing=:none)
    inds = inds[[1,2,3,7,8,9],]
    names!(inds,[:ind_id,:master_id,:skeleton_id,:type,:paper,:date])
#   remove present-day samples
    inds = inds[inds[:date] .!= "0",:]
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
    anno_file = "../../data/ancient_dna/reich_compilation/v37/v37.2.1240K.clean4.anno"
    source_file = "data/aDNA_sources.csv"
    anno = parseAnnoFile(anno_file)
    source = CSV.read(source_file)
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
                 reich_id = join(problems,source,on=:sk_short => :sample,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   joins where my id matched a skeleton id, split by a comma, different order than previous
    for i in 1:nrow(problems)
        problems[i,:sk_short] = split(problems[i,:skeleton_id],", ")[end]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:sample],
                  reich_id = join(problems,source,on=:sk_short => :sample,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   ind_id was a modified version of skeleton id: "_" different position
    for i in 1:nrow(problems)
        problems[i,:sk_short] = split(problems[i,:skeleton_id],"_")[end]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :sk_short,kind=:inner,makeunique=true)[:sample],
                  reich_id = join(problems,source,on=:sk_short=>:sample,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   sample id == :ind_id
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :ind_id,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(problems,source,on=:ind_id => :sample,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   ind_id was a modified version of sample id: "."
    problems[:ind_short] = ""
    for i in 1:nrow(problems)
        problems[i,:ind_short] = split(problems[i,:ind_id],".")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :ind_short,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(problems,source,on=:ind_short=>:sample,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   ind_id was a modified version of sample id: "_"
    for i in 1:nrow(problems)
        problems[i,:ind_short] = split(problems[i,:ind_id],"_")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:sample => :ind_short,kind=:inner,makeunique=true)[:sample],
                  reich_id = join(problems,source,on=:ind_short=>:sample,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
#   sample id was a modified version of :master_id
    source[:id_short] = ""
    for i in 1:nrow(source)
        source[i,:id_short] = split(source[i,:sample],"*")[1]
    end
    id_map = vcat(id_map,DataFrames.DataFrame(colbran_id = join(source,problems,on=:id_short => :master_id,kind=:inner,makeunique=true)[:sample],
                 reich_id = join(problems,source,on=:master_id=>:id_short,kind=:inner,makeunique=true)[:master_id]))
    problems = join(problems,id_map,on=:master_id => :reich_id,kind=:anti)
    println("$(nrow(id_map)) in id_map, $(nrow(problems)) in problems")
end

main()
