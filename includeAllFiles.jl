# Libraries
using JLD

# Create splits Function
include("createSplits.jl")

# best Parameter Function
include("bestParamCrossVal.jl")

# resultsToCSV
include("resultsToCSV.jl")

# Coupling Functions
include("Couplings/bop_Pi.jl")
include("Couplings/cbop_Pi.jl")
include("Couplings/onPath_Pi.jl")
include("Couplings/onPathH_Pi.jl")
include("Couplings/onPathSAst_Pi.jl")
include("Couplings/onPathHSAst_Pi.jl")
include("Couplings/nTimesPath_Pi.jl")
include("Couplings/nTimesPathH_Pi.jl")

# Dissimilarity Functions
include("Dissimilarities/exch_PiToD.jl")
include("Dissimilarities/sur_PiToD.jl")
include("Dissimilarities/freeEnergy_D.jl")

#Kernels Functions
include("Kernels/mds_DToK.jl")
include("Kernels/gauss_DToK.jl")
include("Kernels/markovDiffusion_K.jl")
include("Kernels/modularityMatrix_K.jl")
include("Kernels/regularizedCommuteTime_K.jl")
include("Kernels/cosine_PiToK.jl")

#SClassifMethods
include("SClassifMethods/kernelSVM_Classif.jl")
include("SClassifMethods/sumOfSimilarities_Classif.jl")
include("SClassifMethods/cBoPLabelProp_Classif.jl")
include("SClassifMethods/BoPLabelProp_Classif.jl")
include("SClassifMethods/cBopRS_Classif.jl")
include("SClassifMethods/labelProp_Classif.jl")
include("SClassifMethods/kNN_Classif.jl")

# datasets inclusion
dataNames = [   "mytexas_cocite",
"mywashington_cocite",
"mywisconsin_cocite",
"mycornell_cocite",
"myindustry_pr",
"myindustry_yh",
"mycora",
"myimdb",
"news_2cl_1",
"news_2cl_2",
"news_2cl_3",
"news_3cl_1",
"news_3cl_2",
"news_3cl_3",
"news_5cl_1",
"news_5cl_2",
"news_5cl_3"]
