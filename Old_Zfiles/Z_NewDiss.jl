############################################################
############ Tests Definition
############################################################


# Name your experiment
experimentName = "NewDiss"

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

methodsIndexVec = 1:13

##### Param

betaPV = [1e-4,1e-3,1e-2,1e-1,1,10]
alphaV = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
CP = [0.01, 0.1, 1, 10, 100]

#############################################

# all BoP method
method1 = function(K,y,indexTrain,indexTest,C)
  labelPrediction = kernelSVM_Classif(K,y,indexTrain,indexTest;C=C)
  return(labelPrediction)
end

# METHOD 2 Sum of Sim
method2 = function(A,y,indexTrain,indexTest,alpha)
  labelPrediction = sumOfSimilarities_Classif(A,y,indexTrain,indexTest,alpha)
  return(labelPrediction)
end


#############################################

# METHOD 1 FE gauss
preMethod1 = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 2 Sum of Sim
preMethod2 = function(A)
  return(A)
end

#############################################

# METHOD 3 Modularity
preMethod3 = function(A)
  K = modularityMatrix_K(A)
  return(K)
end

#############################################

# METHOD 4 nonH sur gauss
preMethod4 = function(A,betaP)
  Pi = onPath_Pi(A,betaP)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 5 nonH cosine
preMethod5 = function(A,betaP)
  Pi = onPath_Pi(A,betaP)
  K = cosine_PiToK(Pi)
  return(K)
end

#############################################

# METHOD 6 H sur gauss
preMethod6 = function(A,betaP)
  Pi = onPathH_Pi(A,betaP)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 7 H cosine
preMethod7 = function(A,betaP)
  Pi = onPathH_Pi(A,betaP)
  K = cosine_PiToK(Pi)
  return(K)
end

#############################################

# METHOD 8 ntimes nonH cosine
preMethod8 = function(A,betaP)
  Pi = nTimesPath_Pi(A,betaP)
  K = cosine_PiToK(Pi)
  return(K)
end

#############################################

# METHOD 9 ntimes H cosine
preMethod9 = function(A,betaP)
  Pi = nTimesPathH_Pi(A,betaP)
  K = cosine_PiToK(Pi)
  return(K)
end

#############################################

# METHOD 10 Cov
preMethod10 = function(A,betaP)
  K = onPathCov_K(A,betaP)
  return(K)
end

#############################################

# METHOD 11 Cor
preMethod11 = function(A,betaP)
  K = onPathCor_K(A,betaP)
  return(K)
end

#############################################

# METHOD 12 Cov H
preMethod12 = function(A,betaP)
  K = onPathHCov_K(A,betaP)
  return(K)
end

#############################################

# METHOD 13 Cor H
preMethod13 = function(A,betaP)
  K = onPathHCor_K(A,betaP)
  return(K)
end

################################################
###### METHODS AND PARAMETER STORAGE

allPreParameterValuesArray = [betaPV,false,false,betaPV,betaPV,betaPV,betaPV,betaPV,betaPV,betaPV,betaPV,betaPV,betaPV]
allParameterValuesArray = [CP,alphaV,CP,CP,CP,CP,CP,CP,CP,CP,CP,CP,CP]
allMethods = [method1,method2,method1,method1,method1,method1,method1,method1,method1,method1,method1,method1,method1]
allPreMethods = [preMethod1,preMethod2,preMethod3,preMethod4,preMethod5,preMethod6,preMethod7,preMethod8,preMethod9,preMethod10,preMethod11,preMethod12,preMethod13]
