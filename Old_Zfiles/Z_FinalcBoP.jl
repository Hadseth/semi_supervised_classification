############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "FinalcBoP"

# Include datasets names and index over them
dataIndexVec = 8#[1,3,4,5,6,7,9,10,11,12,13,14,15,16,17]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 5

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 1:9

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
preMethod4 = function(A,betaP)
  Pi = cbop_Pi(A,betaP)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 5 cBopH - Sur - Gauss
preMethod5 = function(A,betaP)
  Pi = cbop_Pi(A,betaP;hitting=true)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 6 cBop w - Sur - Gauss
preMethod6 = function(A,betaP)
  w = sum(A,2)[:] + 1e-50
  w = w / sum(w)
  Pi = cbop_Pi(A,betaP;sIn=w,sOut=w)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 7 cBopH w - Sur - Gauss
preMethod7 = function(A,betaP)
  w = sum(A,2)[:] + 1e-50
  w = w / sum(w)
  Pi = cbop_Pi(A,betaP;sIn=w,sOut=w,hitting=true)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 8 cBop invw - Sur - Gauss
preMethod8 = function(A,betaP)
  w = 1 ./ (sum(A,2)[:] + 1e-50)
  w = w / sum(w)
  Pi = cbop_Pi(A,betaP;sIn=w,sOut=w)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

#############################################

# METHOD 9 cBopH invw - Sur - Gauss
preMethod9 = function(A,betaP)
  w = 1 ./ (sum(A,2)[:] + 1e-50)
  w = w / sum(w)
  Pi = cbop_Pi(A,betaP;sIn=w,sOut=w,hitting=true)
  D = sur_PiToD(Pi)
  K = gauss_DToK(D)
  return(K)
end

################################################
###### METHODS AND PARAMETER STORAGE

allPreParameterValuesArray = [betaPV,false,false,betaPV,betaPV,betaPV,betaPV,betaPV,betaPV]
allParameterValuesArray = [CP,alphaV,CP,CP,CP,CP,CP,CP,CP]
allMethods = [method1,method2,method1,method1,method1,method1,method1,method1,method1]
allPreMethods = [preMethod1,preMethod2,preMethod3,preMethod4,preMethod5,preMethod6,preMethod7,preMethod8,preMethod9]
