############################################################
############ sparse eigen K
############################################################

using LIBLINEAR

######

kernelSVM_Classif = function(A,y,indexTrain,indexTest;nDim=5)

  # make sure its sparse
  Asp = sparse(A)

  # size
  n = size(A)[1]
  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # cost matrix, probability matrix
  C = 1 ./ Asp
  C[!isinf(C)] = C[!isinf(C)]/mean(C[!isinf(C)])
  C[isinf(C)] = myMax
  P = diagm(1 ./ sum(Asp,2)[:,1]) * Asp

  # W and Z matrices
  W = exp(-betaP * C) .* P
  Z = inv(I - W)

  K = 0.5 * (K + K')
  sortind = sortperm(real(eigfact(K)[:values]),rev=true)
  VK = real(eigfact(K)[:vectors])
  VK = VK[:,sortind]
  nDim = convert(Int64,floor(nDim))
  VK = VK[:,1:nDim]

  yf = y[indexTrain]
  VKf = VK[indexTrain,:]
  VKf = (diagm(1 ./sqrt(sum(VKf.*VKf,2)[:]))*(VKf))
  model = linear_train(yf,VKf')

  VKf = VK[indexTest,:]
  VKf = (diagm(1 ./sqrt(sum(VKf.*VKf,2)[:]))*(VKf))
  labPrediction = linear_predict(model,VKf')[1]
  return(labPrediction)
end
