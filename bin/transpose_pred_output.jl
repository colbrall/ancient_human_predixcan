# copied from 1kG_analyses.jl
# transposes predicted_expression.txt files (predixcan output) into gene x ind matrix
# this is the format needed for a lot of the population-level analyses already written
# julia 0.6.1

using DataFrames
using Gadfly
using Cairo
using CSV

# reads file into DataFrame
function readDF(f_path::String)
  println("Reading file $f_path as data frame.......")
  if ispath(f_path)
    df = readtable(f_path,separator='\t', normalizenames=false)
  else
    df = "NA"
  end
  return df
end

# my tissue names -> Eric's tissue names
function mapNames()
  return Dict{String,String}(
    "adipose_subcutaneous" => "Adipose-Subcutaneous",
    "adipose_visceral_omentum" => "adipose_visceral_omentum",
    "brain_putamen_basal_ganglia" => "Brain-Putamen-basalganglia",
    "pancreas" => "Pancreas", "breast_mammary_tissue" => "Breast-MammaryTissue",
    "pituitary" => "Pituitary", "adrenal_gland" => "AdrenalGland",
    "cells_ebv_transformed_lymphocytes" => "Cells-EBV-transformedlymphocytes",
    "anterior_cingulate_cortex" => "Brain-Anteriorcingulatecortex-BA24",
    "cells_transformed_fibroblasts" => "Cells-Transformedfibroblasts",
    "skin_nosun_suprapubic" => "Skin-NotSunExposed-Suprapubic",
    "artery_aorta" => "Artery-Aorta", "colon_sigmoid" => "Colon-Sigmoid",
    "skin_sun_lower_leg" => "Skin-SunExposed-Lowerleg",
    "artery_coronary" => "Artery-Coronary", "colon_transverse" => "Colon-Transverse",
    "small_intestine_terminal_ileum" => "SmallIntestine-TerminalIleum",
    "artery_tibial" => "Artery-Tibial", "spleen" => "Spleen",
    "esophagus_gastroesophageal_junction" => "Esophagus-GastroesophagealJunction",
    "brain_caudate_basal_ganglia" => "Brain-Caudate-basalganglia",
    "esophagus_mucosa" => "Esophagus-Mucosa", "stomach" => "Stomach",
    "brain_cerebellar_hemisphere" => "Brain-CerebellarHemisphere",
    "esophagus_muscularis" => "Esophagus-Muscularis",
    "brain_cerebellum" => "Brain-Cerebellum", "nerve_tibial" => "Nerve-Tibial",
    "heart_atrial_appendage" => "Heart-AtrialAppendage",
    "brain_cortex" => "Brain-Cortex", "liver" => "Liver", "lung" => "Lung",
    "brain_frontal_cortex" => "Brain-FrontalCortex-BA9",
    "brain_hippocampus" => "Brain-Hippocampus", "muscle_skeletal" => "Muscle-Skeletal",
    "whole_blood" => "WholeBlood", "brain_hypothalamus" => "Brain-Hypothalamus",
    "brain_nucleus_accumbens_basal_ganglia" => "Brain-Nucleusaccumbens-basalganglia",
    "prostate" => "prostate","ovary" => "Ovary", "vagina" => "vagina", "testis" => "testis",
    "uterus" => "uterus", "thyroid" => "thyroid","left_ventricle" => "left_ventricle")
end

# converts predicted_expression.txt into an Eric-style file
# takes root directory that holds output directories
function transPredExp(path::String)
  const map_dict = mapNames()::Dict{String,String}
  for item in readdir(path)
    if isdir(realpath("$path/$item"))
      const out_f = ""
      try
        out_f = "$(realpath(path))/$(map_dict[item])_elasticNet0_0.5.full"
      catch
        out_f = "$(realpath(path))/$(item)_elasticNet0_0.5.full"
      end
      if ispath("$(out_f).gz")
        println("$out_f already exists!")
        continue
      end
      const prd_exp = readDF("$path$item/predicted_expression.txt")
      if prd_exp == "NA"
        println("No predicted_expression.txt file in $path$item/")
        continue
      end
      delete!(prd_exp, :FID)
      const genes = names(prd_exp)
      prd_exp = unstack(stack(prd_exp,names(prd_exp)[2:end]),:variable,:IID,:value)
      rename!(prd_exp,:variable => :gene)
      CSV.write(out_f, prd_exp, delim='\t')
      run(`gzip $out_f`)
    end
  end
end


function main()
#### Transposing predicted_expression.txt file
    in_files = ARGS[:,1]::Array{String,1}
    transPredExp(in_files[length(in_files)])

end

main()
