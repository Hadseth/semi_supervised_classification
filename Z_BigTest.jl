############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "BigTest"

# Include datasets names and index over them
dataIndexVec = [18,19,20,5,6]#[1,2,3,4,8,9,10,11,12,13,14,15,16,17]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 5

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 1:32

#############################################
#############################################
#############################################

# Id Method
id_SoS = function(A)
  return(A)
end

###############################

# Modularity K -> modularityMatrix_K(A)

###############################

# FE - MDS K
feMds_K = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = gauss_DToK(D)
  return(K)
end

###############################

# FE - g K
feG_K = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = gauss_DToK(D)
  return(K)
end

###############################

# FE_D -> freeEnergy_D(A,betaP)

###############################

# oPCov -> onPathCov_K(A,betaP)

###############################

# oPCov - SDM
opCovSDM_D = function(A,betaP)
  K = onPathCov_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# oPCovH -> onPathHCov_K(A,betaP)

###############################

# oPCovH - SDM
opCovHSDM_D = function(A,betaP)
  K = onPathHCov_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# oPCor -> onPathCor_K(A,betaP)

###############################

# oPCor - SDM
opCorSDM_D = function(A,betaP)
  K = onPathCor_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# opCorH -> onPathHCor_K(A,betaP)

###############################

# opCorH - SDM
opCorHSDM_D = function(A,betaP)
  K = onPathHCor_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# nTCov -> nTimesPathCov_K(A,betaP)

###############################

# nTCov - SDM
ntCovSDM_D = function(A,betaP)
  K = nTimesPathCov_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# nTCovH -> nTimesPathHCov_K(A,betaP)

###############################

# nTCovH - SDM
ntCovHSDM_D = function(A,betaP)
  K = nTimesPathHCov_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# nTCor -> nTimesPathCor_K(A,betaP)

###############################

# nTCor - SDM
ntCorSDM_D = function(A,betaP)
  K = nTimesPathCor_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

###############################

# nTCorH -> nTimesPathHCor_K(A,betaP)

###############################

# nTCor - SDM
ntCorHSDM_D = function(A,betaP)
  K = nTimesPathHCor_K(A,betaP)
  D = invMDS_KToD(K)
  return(D)
end

#############################################
#############################################
#############################################

# classicalSoS - > sumOfSimilarities_Classif(A,y,indexTrain,indexTest,alpha)

# kSVM
kSVM = function(K,y,indexTrain,indexTest,C)
  labelPrediction = kernelSVM_Classif(K,y,indexTrain,indexTest;C=C)
  return(labelPrediction)
end

##########################

# kSoS
kSoS = function(K,y,indexTrain,indexTest)
  labelPrediction = sumOfK_Classif(K,y,indexTrain,indexTest)
  return(labelPrediction)
end

##########################

# dkNN
kNN = function(D,y,indexTrain,indexTest,kn)
  labelPrediction = kNN_D_Classif(D,y,indexTrain,indexTest;k=kn)
  return(labelPrediction)
end


#############################################
#############################################
#############################################
##################### TESTED METHODS ########################
#
# 1.    SoS
# 2.    Q - SVM
# 3.    Q - SoS
# 4.    FE - MDS - SVM
# 5.    FE - MDS - SoS
# 6.    FE - g - SVM
# 7.    FE - g - SoS
# 8.    FE - kNN
# 9.    oPCov - SVM
# 10.   oPCov - SoS
# 11.   oPCov - SDM - kNN
# 12.   oPCovH - SVM
# 13.   oPCovH - SoS
# 14.   oPCovH - SDM - kNN
# 15.   oPCor - SVM
# 16.   oPCor - SoS
# 17.   oPCor - SDM - kNN
# 18.   oPCorH - SVM
# 19.   oPCorH - SoS
# 20.   oPCorH - SDM - kNN
# 21.   nTCov - SVM
# 22.   nTCov - SoS
# 23.   nTCov - SDM - kNN
# 24.   nTCovH - SVM
# 25.   nTCovH - SoS
# 26.   nTCovH - SDM - kNN
# 27.   nTCor - SVM
# 28.   nTCor - SoS
# 29.   nTCor - SDM - kNN
# 30.   nTCorH - SVM
# 31.   nTCorH - SoS
# 32.   nTCorH - SDM - kNN
#
##########################################

allPreMethods = [   id_SoS,
                    modularityMatrix_K,
                    modularityMatrix_K,
                    feMds_K,
                    feMds_K,
                    feG_K,
                    feG_K,
                    freeEnergy_D,
                    onPathCov_K,
                    onPathCov_K,
                    opCovSDM_D,
                    onPathHCov_K,
                    onPathHCov_K,
                    opCovHSDM_D,
                    onPathCor_K,
                    onPathCor_K,
                    opCorSDM_D,
                    onPathHCor_K,
                    onPathHCor_K,
                    opCorHSDM_D,
                    nTimesPathCov_K,
                    nTimesPathCov_K,
                    ntCovSDM_D,
                    nTimesPathHCov_K,
                    nTimesPathHCov_K,
                    ntCovHSDM_D,
                    nTimesPathCor_K,
                    nTimesPathCor_K,
                    ntCorSDM_D,
                    nTimesPathHCor_K,
                    nTimesPathHCor_K,
                    ntCorHSDM_D]

#

betaPV = [1e-3,1e-2,1e-1,1,10,15]

#

allPreParameterValuesArray = [  false,
                                false,
                                false,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV,
                                betaPV]

allMethods = [                  sumOfSimilarities_Classif,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN,
                                kSVM,
                                kSoS,
                                kNN]

#########

alphaV = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
CP = [0.01, 0.1, 1, 10, 100]
kVec = [1,3,7,13,21]

#########

allParameterValuesArray = [     alphaV,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec,
                                CP,
                                false,
                                kVec]
