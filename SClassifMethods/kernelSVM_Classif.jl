############################################################
############ Kernelized Support Vector machine Kernel
############################################################

using LIBLINEAR

######

kernelSVM_Classif = function(K,y,indexTrain,indexTest;C=1,nDim=5)

  Ks = Symmetric(0.5 * (K + K'))
  n = size(Ks)[1]
  nDim = convert(Int64,floor(nDim))
  decomp = eigfact(Ks,(n-nDim+1):n)

  VK = decomp[:vectors]

  yf = y[indexTrain]
  VKf = VK[indexTrain,:]
  VKf = (diagm(1 ./sqrt.(sum(VKf.*VKf,2)[:]))*(VKf))
  model = linear_train(yf,VKf';eps=1e-10,C=C)

  VKf = VK[indexTest,:]
  VKf = (diagm(1 ./sqrt.(sum(VKf.*VKf,2)[:]))*(VKf))
  labPrediction = linear_predict(model,VKf')[1]

  return(labPrediction)
end
