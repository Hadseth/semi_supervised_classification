############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "NewBigTest"

# Include datasets names and index over them
dataIndexVec = [1,2,3,4,8,9,10,11,12,13,14,15,16,17,18,19,20]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 5

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 1:35

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
# 8.    oPCof - SVM
# 9.    oPCof - SoS
# 10.   oPCofH - SVM
# 11.   oPCofH - SoS
# 12.   oPCos - SVM
# 13.   oPCos - SoS
# 14.   oPCosH - SVM
# 15.   oPCosH - SoS
# 16.   oPCov - SVM
# 17.   oPCov - SoS
# 18.   oPCovH - SVM
# 19.   oPCovH - SoS
# 20.   oPCor - SVM
# 21.   oPCor - SoS
# 22.   oPCorH - SVM
# 23.   oPCorH - SoS
# 24.   nTCos - SVM
# 25.   nTCos - SoS
# 26.   nTCosH - SVM
# 27.   nTCosH - SoS
# 28.   nTCov - SVM
# 29.   nTCov - SoS
# 30.   nTCovH - SVM
# 31.   nTCovH - SoS
# 32.   nTCor - SVM
# 33.   nTCor - SoS
# 34.   nTCorH - SVM
# 35.   nTCorH - SoS
#
##########################################
#############################################
#############################################
#############################################

# Id Method
id_SoS = function(A)
  return(A)
end

###############################

# FE - MDS K
feMds_K = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = mds_DToK(D)
  return(K)
end

###############################

# FE - g K
feG_K = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = gauss_DToK(D)
  return(K)
end

####################################################
####################################################

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

####################################################
####################################################

allPreMethods = [   id_SoS,
                    modularityMatrix_K,
                    modularityMatrix_K,
                    feMds_K,
                    feMds_K,
                    feG_K,
                    feG_K,
                    onPathCof_K,
                    onPathCof_K,
                    onPathHCof_K,
                    onPathHCof_K,
                    onPathCos_K,
                    onPathCos_K,
                    onPathHCos_K,
                    onPathHCos_K,
                    onPathCov_K,
                    onPathCov_K,
                    onPathHCov_K,
                    onPathHCov_K,
                    onPathCor_K,
                    onPathCor_K,
                    onPathHCor_K,
                    onPathHCor_K,
                    nTimesPathCos_K,
                    nTimesPathCos_K,
                    nTimesPathHCos_K,
                    nTimesPathHCos_K,
                    nTimesPathCov_K,
                    nTimesPathCov_K,
                    nTimesPathHCov_K,
                    nTimesPathHCov_K,
                    nTimesPathCor_K,
                    nTimesPathCor_K,
                    nTimesPathHCor_K,
                    nTimesPathHCor_K]

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
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS,
                                kSVM,
                                kSoS]

#########

alphaV = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
CP = [0.01, 0.1, 1, 10, 100]

#########

allParameterValuesArray = [     alphaV,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false,
                                CP,
                                false]
