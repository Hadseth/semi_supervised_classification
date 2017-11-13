############################################################
############ boPLabelProp_Classif
############################################################

using StatsBase

sortsortperm = function(v)
    return sortperm(sortperm(v))
end

kNN_D_Classif = function(D,y,indexTrain,indexTest;k=5)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(D)[1]

  # number of class
  nbc = maximum(y)

  # train set
  redD = D[indexTest,indexTrain]
  maxD = maximum(redD)
  yTrain = y[indexTrain]

  classResult = []

  # "naive" result
  rankD = mapslices(sortsortperm,redD,2)
  knnMat = (rankD .<= k)
  criteriaD = knnMat*maxD - (redD .* knnMat)

  critSumByC = zeros(size(redD)[1],Int64(nbc))
  for cla in 1:Int64(nbc)
      critSumByC[:,cla] = sum(criteriaD[:,yTrain .==cla],2)
  end

  knnRes = mapslices(indmax,critSumByC,2)

  return(knnRes)
end
