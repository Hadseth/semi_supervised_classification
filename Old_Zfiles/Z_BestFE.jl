############################################################
############ Tests Definition
############################################################

# Name your experiment
experimentName = "BestFE"

# Include datasets names and index over them
dataIndexVec = 1#[1,2,3,4,8,9,10,11,12,13,14,15,16,17]

# Number of external/internel folds
nExtFolds = 10
nIntFolds = 5

# Number of Classification
nClassif = 1

#####################################################
#####################################################
###### METHODS AND PARAMETER DEFINITIONS

methodsIndexVec = 1:6

##### Param

betaPV = [1e-3,1e-2,1e-1,1,10,15]
CP = [0.01, 0.1, 1, 10, 100]

#############################################

#############################################
#############################################

# all BoP method
method1 = function(K,y,indexTrain,indexTest,C)
  labelPrediction = kernelSVM_Classif(K,y,indexTrain,indexTest;C=C)
  return(labelPrediction)
end

# all BoP method
method2 = function(K,y,indexTrain,indexTest)
  labelPrediction = sumOfK_classif(K,y,indexTrain,indexTest)
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

# METHOD 2 FE MDS
preMethod2 = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K = mds_DToK(D)
  return(K)
end

#############################################

# METHOD 3 FE gauss deg_norm
preMethod3 = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K_p = gauss_DToK(D)
  Deg = diagm(1 ./ sum(A,2)[:])
  K = Deg*K_p
  return(K)
end

#############################################

# METHOD 4 FE MDS deg_norm
preMethod4 = function(A,betaP)
  D = freeEnergy_D(A,betaP)
  K_p = mds_DToK(D)
  Deg = diagm(1 ./ sum(A,2)[:])
  K = Deg*K_p
  return(K)
end

################################################
###### METHODS AND PARAMETER STORAGE

allPreParameterValuesArray = [betaPV,betaPV,betaPV,betaPV,betaPV,betaPV]
allParameterValuesArray = [CP,CP,false,false,false,false]
allMethods = [method1,method1,method2,method2,method2,method2]
allPreMethods = [preMethod1,preMethod2,preMethod1,preMethod2,preMethod3,preMethod4]
