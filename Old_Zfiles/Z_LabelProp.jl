############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "QuickComp"

# Include datasets names and index over them
dataIndexVec = [1,2,3,4,9,10,11,12,13,14,15,16,17]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 1

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 1:5

##### Param

betaPV = [1e-3,1e-2,1e-1,1,10,15]
alphaV = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]

#############################################

#############################################

# all BoP method
method1 = function(K,y,indexTrain,indexTest)
  labelPrediction = kernelSVM_Classif(K,y,indexTrain,indexTest)
  return(labelPrediction)
end

######

# METHOD 2 Sum of Sim
method2 = function(A,y,indexTrain,indexTest,alpha)
  labelPrediction = sumOfSimilarities_Classif(A,y,indexTrain,indexTest,alpha)
  return(labelPrediction)
end

######

# METHOD 3 Sum of Sim
method3 = function(A,y,indexTrain,indexTest,betaP)
  labelPrediction = sumOfFE_Classif(A,y,indexTrain,indexTest,betaP)
  return(labelPrediction)
end

######

# METHOD 3 Sum of Sim
method4 = function(A,y,indexTrain,indexTest,betaP)
  labelPrediction = sumOfFE2_Classif(A,y,indexTrain,indexTest,betaP)
  return(labelPrediction)
end

######

# METHOD 3 Sum of Sim
method5 = function(A,y,indexTrain,indexTest,betaP)
  labelPrediction = sumOfFE3_Classif(A,y,indexTrain,indexTest,betaP)
  return(labelPrediction)
end

#############################################

# METHOD 1 OnPathPi old
preMethod1 = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = gauss_DToK(D)
  return(K)
end

####################################

# METHOD 2 Sum of FE
preMethod2 = function(A)
  return(A)
end

################################################
###### METHODS AND PARAMETER STORAGE

allPreParameterValuesArray = [betaPV,false,false,false,false]
allParameterValuesArray = [false,alphaV,betaPV,betaPV,betaPV]
allMethods = [method1,method2,method3,method4,method5]
allPreMethods = [preMethod1,preMethod2,preMethod2,preMethod2,preMethod2]
