############################################################
############ boPLabelProp_Classif
############################################################

cBopRS_Classif = function(A,y,indexTrain,indexTest,betaP;weights=repeat([1/(size(A)[1])],inner=size(A)[1]),hitting = false,persist=1,convThres=1e-10)

  # myMax, myMin
  myMin   = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(A)[1]

  # number of class
  nbc = Int64(maximum(y))

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)

  # cost matrix, probability matrix
  C = 1 ./ A
  C[!isinf(C)] = C[!isinf(C)]/mean(C[!isinf(C)])
  C[isinf(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:,1]) * A

  answMat = zeros(nbc,sum(indexTest))
  for class in 1:nbc
    if sum((indexTrain)&(y .== class)) != 0
      # supply and demand
      sIn = zeros(n)
      sIn[(indexTrain)&(y .== class)] = weights[(indexTrain)&(y .== class)]
      sIn = sIn / sum(sIn)
      sOut = sIn

      #if hitting
      #  W = exp(-betaP * C) .* P
      #  Z = inv(I - W)
      #  divisor = diag(Z)
      #  numerator = divisor
      #  outVal = sOut
      #else
      Q = I - P'
      pinvQ = pinv(Q)
      q = (I - pinvQ * Q)
      q = q[:,1]

      nRef0 = pinvQ * (sIn - P' * sOut)
      gammaP = maximum((sOut - nRef0) ./ q) + persist
      nRef = nRef0 + gammaP * q
      alpha = sOut ./ nRef
      Pmod = (I - diagm(alpha)) * P

      # matrices W and Z, vectors divisor and out_val
      W = exp(-betaP * C) .* Pmod
      Z = inv(I - W)
      divisor = ones(n)
      numerator = nRef
      outVal = alpha
      #end

      # Iterations
      mIn = e
      mOut = e
      convergence = false
      FE = myMax
      while !convergence
        FEprev = FE
        mIn = e ./ (Z * (mOut .* (outVal ./ divisor)))
        mOut = numerator ./ (Z' * (mIn .* sIn))
        FE = sum(-log(mIn)/betaP .* sIn) + sum(-log(mOut)/betaP .* sOut)
        maxDiff = abs(FE - FEprev)
        if maxDiff < convThres
          convergence = true
        end
      end
      answMat[class,:] = (nRef .* ((e ./ (mIn .* mOut)) - alpha))[indexTest]
    end
  end
  yp = findmax(answMat,1)[2]
  labelPrediction = ((yp-1)%nbc + 1)[:]
  return(labelPrediction)
end
