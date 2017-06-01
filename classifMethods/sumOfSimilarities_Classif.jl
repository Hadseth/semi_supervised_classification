############################################################
############ Sum of Similarities Classifier
############################################################

sumOfSimilarities_Classif = function(A,y,indexTrain,indexTest,alpha;convThres=1e-5)

  n = length(y)

  D0 = diagm(1 ./ sum(A,2)[:])
  # transition probability matrix
  P = D0 * A
  # number of classes
  nbc = maximum(y)

  # only the training samples are labeled the others are set to a dummy class 0
  yt  = zeros(n,1)
  yt[indexTrain] = y[indexTrain]

  Y1  = Array{Float64}(n,0)
  for c = 1:nbc
    Y1 = [Y1 (yt.==c)]
  end

  # predisted classes initialized to class memberships
  Yp = Y1
  convergence = false
  while !convergence
    YpPrev = Yp
    # computing the sum over forests from the randomized commute-time kernel
	  Yp = alpha*P'*Yp + (1-alpha)*Y1
    if vecnorm(YpPrev - Yp) < convThres
      convergence = true
    end
  end
  Yp = D0 * Yp

  # assign to class with maximal similarity
  yp = findmax(Yp',1)[2]
  yp = (yp-1)%nbc + 1
  yp = yp'

  return(yp[indexTest])
end
