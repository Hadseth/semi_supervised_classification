############################################################
############ Contrained BoP Coupling
############################################################

cbopSp_Pi = function(A,betaP;sIn=repeat([1/(size(A)[1])],inner=size(A)[1]),sOut=repeat([1/(size(A)[1])],inner=size(A)[1]),persist=0.1,convThres=1e-10)

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
  P = sparse(diagm(1 ./ sum(A,2)[:,1]) * A)

  # supply and demand
  sIn = sIn / sum(sIn)
  sOut = sOut / sum(sOut)

  W = exp.(-betaP * C) .* P
  Z = inv(I - W)
  Zh = Z * diagm(1 ./ diag(Z))
  invZh = inv(Zh)
  invZh = diagm(diag(Z))*(I - W)

  # Iterations
  imIn = e
  imOut = e
  FE = myMax
  convergence = false

  while !convergence
    imInPrev = imIn

    imOut = sOut ./ (invZh*imIn)
    imIn = sIn ./ (invZh'*imOut)
    maxDiff = maximum(abs(imInPrev - imIn))
    if maxDiff < convThres
      convergence = true
    end
  end
  Pi2 = diagm(sIn ./ imIn) * Zh * diagm(sOut ./ imOut)

  # Iterations
  mIn = e
  mOut = e
  convergence = false
  FE = myMax
  while !convergence
    FEPrev = FE
    mIn = 1 ./ (Zh * (sOut .* mOut))
    mOut = 1 ./ (Zh' * (sIn .* mIn))
    FE = sum(-log.(mIn) .* sIn) + sum(-log.(mOut) .* sOut)
    maxDiff = abs(FEPrev - FE)
    if maxDiff < convThres
      convergence = true
    end
  end
  Pi = diagm(sIn .* mIn) * Zh * diagm(sOut .* mOut)

  return(Pi)
end
