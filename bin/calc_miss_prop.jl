# calc_miss_prop.jl
#
#
# julia1.1

using DataFrames
using CSV
using ArgParse
using GZip

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--outf","-o"
            help = "output file"
            arg_type = String
            default = "miss_prop.txt"

        "--addMiss","-a"
            help = "number of additional variants to consider 'missing' that aren't in the dosage file at all"
            arg_type = Int64
            default = 0

        "--pop_dir","-p"
            help = "directory of dosage files to count missingness for"
            arg_type = String

        "--samples","-s"
            help = "sample.txt file to go with pop_dir"
            arg_type = String

        "--postfix","-f"
            help = "ending for dosage files. Default is .dos.gz"
            arg_type = String
            default = ".dos.gz"
    end
    return parse_args(s)
end

function calcMiss(pop_dir::String, samples::String, miss::Int64, outf::String,postfix::String)
    df = CSV.read(samples; delim = "\t", header=[:ind_id, :fam_id])
    deletecols!(df,:fam_id)
    num_miss = fill(0,nrow(df))
    total_snps = miss #zero if I shouldn't be counting any snps not in these dosage files.
    for item in readdir(pop_dir)
        if !endswith(item, postfix) continue end
        println(item)
        GZip.open("$pop_dir$item","r") do inf
            for line in eachline(inf)
                l = split(chomp(line),"\t")[7:end]
                miss_inds = findall(i->(!isequal(i,"0") && !isequal(i,"1") && !isequal(i,"2")),l)
                num_miss[miss_inds] .+= 1
                total_snps += 1
            end
        end
    end
    println("Total SNPs = $total_snps")
    df[:prop_miss] = num_miss/total_snps
    CSV.write(outf,df;delim="\t")
end

function main()
    args = parse_commandline()
    calcMiss(args["pop_dir"],args["samples"],args["addMiss"],args["outf"],args["postfix"])
end

main()
