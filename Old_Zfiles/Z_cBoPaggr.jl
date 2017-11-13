############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "cBoPaggr"

# Include datasets names and index over them
dataIndexVec = [1,3,4,8,9,10,11,12,13,14,15,16,17]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 5

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 4:5#1:5

##### Param

betaPV = [1e-4,1e-3,1e-2,1e-1,1,10]
powV = [-1,-0.5,0,0.5,1]
paramsV = [ [par1,par2]  for par1 in betaPV, par2 in powV]
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

# METHOD 1 FE - GAUSS
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

# METHOD 4 cBop - Sur - Gauss
preMethod4 = function(A,betaP,pow)
  w = sum(A,2)[:] + 1e-50
  wp = (w).^pow
  Pi = cbop_Pi(A,betaP;sIn=wp,sOut=wp)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 5 cBop - Sur - Gauss
preMethod5 = function(A,betaP,pow)
  w = sum(A,2)[:] + 1e-50
  wp = (w).^pow
  Pi = cbop_Pi(A,betaP;sIn=wp,sOut=wp,hitting=true)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

################################################
###### METHODS AND PARAMETER STORAGE

allPreParameterValuesArray = [betaPV,false,false,paramsV,paramsV]
allParameterValuesArray = [CP,alphaV,CP,CP,CP]
allMethods = [method1,method2,method1,method1,method1]
allPreMethods = [preMethod1,preMethod2,preMethod3,preMethod4,preMethod5]
