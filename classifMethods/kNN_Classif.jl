############################################################
############ boPLabelProp_Classif
############################################################

using StatsBase

kNN_Classif = function(K,y,indexTrain,indexTest;k=5,wei = false)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(K)[1]

  # number of class
  nbc = maximum(y)

  # train set
  redK = K[indexTest,indexTrain]
  yTrain = y[indexTrain]

  classResult = []
  weightsMatrix = copy(redK)
  if !wei
    weightsMatrix[:,:] = 1
  end

  for rind in 1:sum(indexTest)
    rowK = redK[rind,:]
    maxK = sortperm(sortperm(rowK)) .> sum(indexTrain) - k
    class = indmax(counts(Array{Int64}(yTrain[maxK]),Int64(nbc),weights(weightsMatrix[rind,maxK])))
    push!(classResult,class)
  end

  return(classResult)
end
