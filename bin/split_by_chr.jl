#split_by_chr.jl
# @author Laura Colbran
#
# julia 1.1

using ArgParse
using GZip

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--vcf"
            help = "if input is vcf file"
            action = :store_true
        "--files","-f"
            nargs= '*'
            help = "files(s) to split by chr"
            arg_type = String
    end
    return parse_args(s)
end

function splitVCF(files::Array{String,1})
    for f in files
        out_fns = ["$(split(f,".")[1])_chr$(i).vcf" for i in 1:22]
        outfs = [open(g, "w") for g in out_fns]
        GZip.open(f,"r") do file
            for line in eachline(file)
                #println(line)
                if startswith(line,"##contig")
                    chr = tryparse(Int64,split(split(line,"=")[3],",")[1])
                    if typeof(chr) == Nothing continue end
                    write(outfs[chr],"$(line)\n")
                elseif startswith(line,"#")
                    for i in 1:22
                        write(outfs[i],"$(line)\n")
                    end
                else
                    chr = tryparse(Int64,split(line, "\t")[1])
                    if typeof(chr) == Nothing continue end
                    write(outfs[chr],"$(line)\n")
                end
            end
        end
        for i in 1:22
            close(outfs[i])
        end
    end
end

function main()
    parsed_args = parse_commandline()
    if parsed_args["vcf"]
        splitVCF(parsed_args["files"])
    end
end

main()
