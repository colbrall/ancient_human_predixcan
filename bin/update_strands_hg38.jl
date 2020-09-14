# update_strands_hg38.jl
#
# updates the ref/alt alleles that had strand flips between hg19 and hg38 by referencing the 1kG hg38 VCF files.
# Also pulls pop-level AF info, and (corrected) Ancestral allele calls from the 1kG VCFs,
# and sorts out any duplicate lines that in the input that were due to ref/alt issues
# assumes you're giving it chr-specific files (can be all or a specific one).

using ArgParse
using GZip

const CHRS = ["22","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","X","Y"]
const COMPLEMENT = Dict{String,String}("A"=>"T","C"=>"G","G"=>"C","T"=>"A","a"=>"t","c"=>"g","g"=>"c","t"=>"a","."=>".","N"=>"N")
const INFO_VARS = ["AF","EAS_AF","EUR_AF","AMR_AF","AFR_AF","SAS_AF","AA"] # ones to pull from 1kG
const INFO_INDS = [2,6,7,8,9,10,11] # indices in 1kG info column for variables specified in INFO_VARS
const FLIP_FLAG = "STRAND_FLIP" # flag in 1kG to indicate strand was flipped

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--genome_path","-g"
            help = "path to files containing genome I'm wanting to fix. assumes gzipped vcf files split by chromosome. Put * in place of chr number (include quotes around it!), and script will fill it in"
            arg_type = String

        "--ref_path","-r"
            help = "directory with 1kG vcfs. structure same as genome_path."
            arg_type = String

        "--out_dir", "-o"
            help = "directory path to write output files"
            arg_type = String
            default = "./"
    end
    return parse_args(s)
end

function readGen(path::String)
    dict = Dict{Tuple{Int64,Int64},Array{String,1}}() #(CHR,POS) => [ID,REF,ALT,QUAL|FILT,INFO,FORM|GENOS]
    header = String[]
    GZip.open(path) do f
        open("duplicates.txt","w") do dup
            for line in eachline(f)
                if startswith(line,"#")
                    append!(header,[chomp(line)])
                    continue
                end
                l = split(chomp(line),"\t")
                snp = (parse(Int64,l[1]),parse(Int64,l[2]))
                # handle duplicates
                if haskey(dict, snp)
                    if l[5] == dict[snp][2] && dict[snp][3] == "."
                        dict[snp][2] = l[4]
                        dict[snp][3] = l[5]
                        gts = split(dict[snp][5],"\t")
                        for i in 1:length(gts)
                            if gts[i] == "0/0"
                                gts[i] = "1/1"
                            elseif gts[i] == "./." #fill in missing genotypes
                                gts[i] = l[8+i]
                            end
                        end
                        dict[snp][5] = join(gts,"\t")
                    elseif l[5] == "." && dict[snp][3] == "."
                        dict[snp][3] = l[5]
                        gts = split(dict[snp][5],"\t")
                        for i in 1:length(gts)
                            if gts[i] == "./."
                                if l[8+i] == "0/0"
                                    gts[i] = "1/1"
                                else
                                    gts[i] = l[8+i]
                                end
                            end
                        end
                        dict[snp][5] = join(gts,"\t")
                    else
                        write(dup,line)
                        write(dup,"$(dict[join(l[1:2],"\t")])\n")
                    end
                end
                dict[snp] = [l[3],l[4],l[5],join(l[6:7],"\t"),l[8],join(l[9:end],"\t")]
            end
        end
    end
    return dict,header
end

function updateStrand(gen_path::String,ref_path::String,out_dir::String)
    by_chr = true
    if !occursin("*",gen_path)
        by_chr = false
    end
    for chr in CHRS
        inf = replace(gen_path,"*"=>chr)
        if !ispath(inf)
            println("chr$chr not in genome")
            continue
        end
        gen_dict,header = readGen(inf) #(CHR,POS) => [ID,REF,ALT,QUAL|FILT,INFO,FORM|GENOS]
        out_path = "$(out_dir)$(split(basename(inf),".gz")[1])"
        head_writ = false
        open(out_path,"w") do outf
            for line in header[1:end-1]
                write(outf,"$line\n")
            end
            GZip.open(replace(ref_path,"*"=>chr)) do f
                for line in eachline(f)
                    if startswith(line,"##INFO")
                        var = "$(split(split(chomp(line),",")[1],"=")[end])"
                        if in(var,INFO_VARS)
                            write(outf,"$line\n")
                        end
                        continue
                    elseif startswith(line,"#")
                        continue
                    end
                    if !head_writ
                        if length(header) > 0 write(outf,"$(header[end])\n") end
                        head_writ = true
                    end
                    snp = (parse(Int64,split(line,"\t")[1]),parse(Int64,split(line,"\t")[2]))
                    if !haskey(gen_dict,snp) continue end
                    ref_info = split(split(line,"\t")[8],";")
                    # match ref/alt so any 1kG allele frequencies will be right
                    if split(line,"\t")[4] != gen_dict[snp][2]
                        ref = gen_dict[snp][3]
                        alt = gen_dict[snp][2]
                        tmp = split(gen_dict[snp][6],"\t")
                        for i in 2:length(tmp)
                            if tmp[i] == "0/0"
                                tmp[i] = "1/1"
                            elseif tmp[i] == "0/1"
                                tmp[i] = "1/0"
                            elseif tmp[i] == "1/0"
                                tmp[i] = "0/1"
                            elseif tmp[i] == "1/1"
                                tmp[i] = "0/0"
                            end
                        end
                        gen_dict[snp][6] = join(tmp,"\t")
                    else
                        ref = gen_dict[snp][2]
                        alt = gen_dict[snp][3]
                    end
                    # for each snp, add info of interest to dict info
                    out_info = append!(split(gen_dict[snp][5],";"),ref_info[INFO_INDS])
                    # fix STRAND_FLIP and ancestral call if necessary
                    if ref_info[end] == FLIP_FLAG
                        ref = COMPLEMENT[ref]
                        alt = COMPLEMENT[alt]
                        if in("AA",INFO_VARS)
                            ind = findfirst([startswith(x,"AA=") for x in out_info])
                            #println(out_info[ind])
                            if typeof(ind) != Nothing
                                aa = split(split(out_info[ind],"=")[2],"|")
                                for i in 1:length(aa)
                                    try
                                        aa[i] = COMPLEMENT[aa[i]]
                                    catch
                                        continue
                                    end
                                end
                                out_info[ind] = "AA=$(join(aa,'|'))"
                            end
                        end
                    end
                    gen_dict[snp][2] = ref  #CHR|POS => [ID,REF,ALT,QUAL|FILT,INFO,FORM|GENOS]
                    gen_dict[snp][3] = alt
                    gen_dict[snp][5] = join(out_info,";")
                end
            end
            for key in sort(collect(keys(gen_dict)))
                println(key)
                write(outf,"$(key[1])\t$(key[2])\t$(gen_dict[key][1])\t$(gen_dict[key][2])\t$(gen_dict[key][3])\t$(gen_dict[key][4])\t$(gen_dict[key][5])\t$(gen_dict[key][6])\n")
            end
        end
        run(`bgzip -f $out_path`)
        if !by_chr break end
    end
end

function main()
    parsed_args = parseCommandLine()
    updateStrand(parsed_args["genome_path"],parsed_args["ref_path"],parsed_args["out_dir"])
end

main()
