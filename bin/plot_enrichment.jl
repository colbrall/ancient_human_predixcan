# plot_enrichment.jl
# @author Laura Colbran
#
# plots enrichment results from WebGestalt output
#
# julia1.1

using ArgParse
using DataFrames
using CSV
using Seaborn

#set plotting defaults
# default(color=:black,leg=false,grid=false,fontfamily="arial")
Seaborn.set(style="white", palette="muted")
set_style(Dict("font.family" =>["DejaVu Sans"]))

function parseCommandLine()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--enr_file","-f"
            help = "path to WebGestalt output file"
            arg_type = String
            required = true
    end
    return parse_args(s)
end

function main()
    parsed_args = parseCommandLine()
    enrichments = CSV.read(parsed_args["enr_file"])
    enrichments[:log] = [-log(10,x) for x in enrichments[Symbol("P Value")]]
    # println(first(enrichments,6))

    e_plot = Seaborn.scatter(enrichments[:Ratio],enrichments[:log],sizes=enrichments[:Size],
                color=:black)
    Seaborn.xlabel("O/E Enrichment Ratio")
    Seaborn.ylabel("-log10(P Value)")
    # Seaborn.despine()
    Seaborn.savefig("$(split(basename(parsed_args["enr_file"]),".")[1])_scatter.pdf")
end

main()
