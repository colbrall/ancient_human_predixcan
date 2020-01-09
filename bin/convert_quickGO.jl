# for converting from quickGO gene output to Ensembl IDs and locations
#
# likely to be a combo of generalization and specific ID matching
#

using ArgParse

#IDs manually mapped using ensembl.org
const ID_MAP = Dict{String,String}("A0A384N679"=> "CPE","IHPK1"=>"IP6K1","A0A140VK19"=>"GATM","A0A140VK13"=> "CPT2","DGIC"=>"PLAC8","TNLG2A"=>"TNFSF18",
                    "ZNF278"=>"PATZ1","HEL-S-2a"=>"PRDX2","CREA7-4"=>"NKX2-4","ARHH"=>"RHOH","AMICA1"=>"JAML","EBI2"=>"GPR183","MNAB"=>"RC3H2",
                    "ATL1-alpha"=>"BCL11B","H0Y858"=>"AC097637.1","STING"=>"TMEM173","ILRUN"=>"C6orf106","LGP2"=>"DHX58","CGAS"=>"MB21D1","MITA"=>"TMEM173",
                    "HEL-S-19"=>"GALK1","NUP42"=>"NUPL2","PRXL2C"=>"AAED1","HEL-S-87p"=>"ALDOA","HEL-S-49"=>"TPI1","HEL-S-278"=>"GAPDHS","HEL-S-162eP"=>"GAPDH",
                    "PKM2"=>"PKM","A0A140VJR3"=>"PGK2","HEL-S-30"=>"PKM")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--parsed_file","-p"
            help = "predixcan-parsed file with gene lists"
            arg_type = String
        "--go_file","-g"
            nargs='*'
            help = "quickGO .tsv file(s) to convert"
            arg_type = String
    end
    return parse_args(s)
end

function geneDict(path::String)
    genes = Dict{SubString{String},Array{SubString{String},1}}()
    open(path) do f
        for line in eachline(f)
            l = split(chomp(line),"\t")
            genes[l[3]] = l[[1,4,5,2]]
        end
    end
    return genes
end

function collapseSet(path::String)
    genes = Dict{SubString{String},Set{SubString{String}}}()
    open(path) do f
        for line in eachline(f)
            l = split(chomp(line),"\t")
            if in(l[3],keys(genes))
                push!(genes[l[3]],l[5])
            else
                genes[l[3]] = Set([l[5]])
            end
        end
    end
    return genes
end

function convertGO(go_paths::Array{String,1},ref_path::String)
    ref_genes = geneDict(ref_path)
    for file in go_paths
        out_p = splitext(file)[1]
        targets = collapseSet(file)
        open("$(out_p)_ids.bed","w") do success
            write(success,"#chr\tstart\tend\tensembl_id\tname\tGO_terms\n")
            open("$(out_p)_nonconvert.txt","w") do fail
                for gene in keys(targets)
                    if in(gene,keys(ref_genes))
                        loc = join(ref_genes[gene],"\t")
                        terms = join(targets[gene],";")
                        write(success,"$loc\t$gene\t$terms\n")
                    elseif in(gene,keys(ID_MAP))
                        loc = join(ref_genes[ID_MAP[gene]],"\t")
                        terms = join(targets[gene],";")
                        write(success,"$loc\t$gene\t$terms\n")
                    else
                        write(fail,"$gene\n")
                    end
                end
            end
        end
    end
end

function main()
    parsed_args = parse_commandline()
    convertGO(parsed_args["go_file"],parsed_args["parsed_file"])
end

main()
