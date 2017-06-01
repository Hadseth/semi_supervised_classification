############################################################
############ Kernelized Support Vector machine Kernel
############################################################

using LIBLINEAR

######

kernelSVM_Classif = function(K,y,indexTrain,indexTest;nDim=5)
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
