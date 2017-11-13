############################################################
############ Contrained BoP Coupling
############################################################

cbop_Pi = function(A,betaP;hitting=false,sIn=repeat([1/(size(A)[1])],inner=size(A)[1]),sOut=repeat([1/(size(A)[1])],inner=size(A)[1]),persist=0.1,convThres=1e-10)

  # myMax, myMin
  myMin = 1e-30
  myMax = typemax(Int64)/1000

  # size
  n = size(A)[1]

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)

  # cost matrix, probability matrix
  C = 1 ./ A
  C[.!isinf.(C)] = C[.!isinf.(C)]/mean(C[.!isinf.(C)])
  C[isinf.(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:,1]) * A

  # supply and demand
  sIn = sIn / sum(sIn)
  sOut = sOut / sum(sOut)

  if hitting
    W = exp.(-betaP * C) .* P
    Z = inv(I - W)
    Z = Z * diagm(1 ./ diag(Z))
    outVal = sOut
    numerator = e
  else
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
    W = exp.(-betaP * C) .* Pmod
    Z = inv(I - W)
    outVal = alpha
    numerator = nRef
  end

  # Iterations
  mIn = e
  mOut = e
  convergence = false
  FE = myMax
  while !convergence
    FEPrev = FE
    mIn = e ./ (Z * (mOut .* outVal))
    mOut = numerator ./ (Z' * (mIn .* sIn))
    FE = sum(-log.(mIn) .* sIn) + sum(-log.(mOut) .* sOut)
    maxDiff = abs(FEPrev - FE)
    if maxDiff < convThres
      convergence = true
    end
  end
  Pi = diagm(mIn .* sIn) * Z * diagm(mOut .* outVal)

  return(Pi)
end
