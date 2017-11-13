############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "QuickComp"

# Include datasets names and index over them
dataIndexVec = 20#[1,2,3,4,8,9,10,11,12,13,14,15,16,17]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 1

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 1:3

##### Param

betaPV = [1e-3,1e-2,1e-1,1,10,15]
alphaV = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
CP = [0.01, 0.1, 1, 10, 100]
kVec = [1,3,7,13,21]

#############################################

#############################################
#############################################

# all BoP method
method1 = function(A,y,indexTrain,indexTest,alpha)
  labelPrediction = sumOfSimilarities_Classif(A,y,indexTrain,indexTest,alpha)
  return(labelPrediction)
end

# all BoP method
method2 = function(K,y,indexTrain,indexTest)
  labelPrediction = kernelSVM_Classif(K,y,indexTrain,indexTest)
  return(labelPrediction)
end

method3 = function(K,y,indexTrain,indexTest)
  labelPrediction = sumOfK_Classif(K,y,indexTrain,indexTest)
  return(labelPrediction)
end

#############################################

# METHOD 1 FE gauss
preMethod1 = function(A)
  return(A)
end

#############################################

# METHOD 1 non-hitting
preMethod2 = function(A, betaP)
  D = freeEnergy_D(A, betaP)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 1 non-hitting
preMethod3 = function(A, betaP)
  D = freeEnergy_D(A, betaP)
  K = mds_DToK(D)
  return(K)
end


################################################
###### METHODS AND PARAMETER STORAGE

allPreParameterValuesArray = [false, betaPV, betaPV]
allParameterValuesArray = [alphaV, false, false]
allMethods = [method1, method2, method3]
allPreMethods = [preMethod1, preMethod2, preMethod3]
