############################################################
############ boPLabelProp_Classif
############################################################

labelProp_Classif = function(K,y,indexTrain,indexTest)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(K)[1]

  # number of class
  nbc = maximum(y)

  # train set
  yt  = zeros(n,1)
  yt[indexTrain] = y[indexTrain]

  # Y all
  Y1  = Array{Float64}(n,0)
  for class = 1:nbc
    colY1 = (yt.==class) - ((yt.!=class) & (yt.!= 0))
    Y1 = [Y1 colY1]
  end

  Yp = K * Y1

  # assign to class with maximal similarity
  yp = findmax(Yp',1)[2]
  yp = (yp-1)%nbc + 1
  yp = yp'

  return(yp[indexTest])
end
