# convert dbSNP bed file to annotation file for PrediXcan (to go into preprocessing)
#
# julia 1.1

using ArgParse
using GZip
using DataFrames
using CSV

# parses command-line arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--bed_file","-b"
            help = "dbSNP bed file to convert"
            arg_type = String
        "--out", "-o"
            help = "output file prefix"
            arg_type = String
    end
    return parse_args(s)
end

function bedToAnno(file_path::String,out_name::String)
    GZip.open(file_path,"r") do inf
        open("$out_name.anno.txt","w") do outf
            write(outf,"#Chr\tPos\tVariantID\tRef\tAlt\trsID")
            for line in eachline(inf)
                l = split(chomp(line),"\t")
                chr = split(l[1],"r")[end]
                alleles = split(l[9],"/")
                alt = ""
                ref = ""
                if length(alleles) > 2   #skip anything that's not biallelic
                    println(line)
                    exit()
                elseif length(alleles) < 2 #fill in alleles for indels
                    if l[7] == "-"
                        alt = split(l[22],",")[2]
                        ref = split(alt,"")[1]
                    else
                        ref = split(l[22],",")[2]
                        alt = split(ref,"")[1]
                    end
                else
                    if alleles[1] == l[7]
                        alt = alleles[2]
                        ref = alleles[1]
                    else
                        alt = alleles[1]
                        ref = alleles[2]
                    end
                    #fill in alleles for other indels
                    if ref == "-"
                        ref = split(alt,"")[1]
                    elseif alt == "-"
                        alt = split(ref,"")[1]
                    end
                end
                write(outf,join([chr,l[3],"$(chr)_$(l[3])_$(ref)_$(alt)_b38",ref,alt,"$(l[4])\n"],"\t"))
            end
        end
    end
end

function main()
    parsed_args = parse_commandline()
    bedToAnno(parsed_args["bed_file"],parsed_args["out"])
end

main()
