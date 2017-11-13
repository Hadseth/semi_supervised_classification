############################################################
############ Sum of Similarities Classifier
############################################################

sumOfFE_Classif = function(A,y,indexTrain,indexTest,betaP;convThres=1e-3)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(A)[1]

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)
  H = I - e*e' / n

  # cost matrix, probability matrix
  C = 1 ./ A
  C[.!isinf.(C)] = C[.!isinf.(C)]/mean(C[.!isinf.(C)])
  C[isinf.(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:,1]) * A

  # W and Z matrices
  W = exp.(-betaP * C) .* P
  alpha = 1 .- sum(W,2)[:]

  # number of classes
  nbc = maximum(y)

  # only the training samples are labeled the others are set to a dummy class 0
  yt  = zeros(n,1)
  yt[indexTrain] = y[indexTrain]

  Y1  = Array{Float64}(n,0)
  for c = 1:nbc
    Y1 = [Y1 (yt.==c)]
  end
  Ynorm = Y1 * diagm(1 ./ (sum(Y1,1)[:]+1e-40))
  sigma = sum(Ynorm,2)[:]
  # predisted classes initialized to class memberships
  Yp = Y1
  convergence = false
  while !convergence
    YpPrev = Yp
    # computing the sum over forests from the randomized commute-time kernel
	  Yp = W'*Yp + ((sigma*(alpha'*Yp)).* Y1)
    if vecnorm(YpPrev - Yp) < convThres
      convergence = true
    end
  end

  # assign to class with maximal similarity
  yp = findmax(Yp',1)[2]
  yp = (yp-1).%nbc + 1
  yp = yp'

  return(yp[indexTest])
end
