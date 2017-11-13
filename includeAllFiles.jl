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

# Dissimilarity Functions
include("Dissimilarities/exch_PiToD.jl")
include("Dissimilarities/sur_PiToD.jl")
include("Dissimilarities/freeEnergy_D.jl")
include("Dissimilarities/approxFreeEnergy_D.jl")
include("Dissimilarities/invMDS_KToD.jl")
include("Dissimilarities/nodeDistrib_D.jl")


#Kernels Functions
include("Kernels/mds_DToK.jl")
include("Kernels/gauss_DToK.jl")
include("Kernels/markovDiffusion_K.jl")
include("Kernels/modularityMatrix_K.jl")
include("Kernels/regularizedCommuteTime_K.jl")
include("Kernels/cosine_PiToK.jl")
include("Kernels/onPathCor_K.jl")
include("Kernels/onPathCov_K.jl")
include("Kernels/onPathCos_K.jl")
include("Kernels/onPathCof_K.jl")
include("Kernels/onPathHCor_K.jl")
include("Kernels/onPathHCov_K.jl")
include("Kernels/onPathHCos_K.jl")
include("Kernels/onPathHCof_K.jl")
include("Kernels/nTimesPathCor_K.jl")
include("Kernels/nTimesPathCov_K.jl")
include("Kernels/nTimesPathCos_K.jl")
include("Kernels/nTimesPathHCor_K.jl")
include("Kernels/nTimesPathHCov_K.jl")
include("Kernels/nTimesPathHCos_K.jl")


#SClassifMethods
include("SClassifMethods/kernelSVM_Classif.jl")
include("SClassifMethods/sumOfSimilarities_Classif.jl")
include("SClassifMethods/cBoPLabelProp_Classif.jl")
include("SClassifMethods/BoPLabelProp_Classif.jl")
include("SClassifMethods/kNN_D_Classif.jl")
include("SClassifMethods/sumOfK_Classif.jl")
include("SClassifMethods/sumOfFE_Classif.jl")

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
"news_5cl_3",
"email",
"blogs",
"polbooks"]
