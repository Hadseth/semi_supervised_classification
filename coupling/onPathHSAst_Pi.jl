############################################################
############ BoP Coupling
############################################################

onPathHSAst_Pi = function(A,betaP)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(A)[1]

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)

  # cost matrix, probability matrix
  C = 1 ./ A
  C[!isinf(C)] = C[!isinf(C)]/mean(C[!isinf(C)])
  C[isinf(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:]) * A

  # W and Z matrices
  W = exp(-betaP * C) .* P
  Z = inv(I - W)
  da = diag(Z)

  Zh = Z ./ (e*da')

  Pi = zeros(n,n)
  for t in 1:n
    Zhittj = (Zh[:,t] * Zh[t,:]')
    ZhWOt = (Zh -  Zhittj) * diagm(1 ./ (1 - diag(Zhittj)))
    ZhWOt[:,t] = 0
    denomZhWOt = (1 - ZhWOt.*ZhWOt')
    Endt = repmat(Zh[:,t]',n,1)
    PlusPi = @parallel (+) for s in 1:n
      ZhsiWOtMat = repmat(ZhWOt[s,:],1,n)
      ZFiSjWOt = ((ZhsiWOtMat - (ZhsiWOtMat.*ZhWOt)') ./ denomZhWOt)
      ZFiSjWOt[diagind(ZFiSjWOt)] = 0

      addPiNotSym = (ZFiSjWOt  .* ZhWOt .* Endt) / Zh[s,t]
      addPi = addPiNotSym + addPiNotSym'
      specialCaseVec = ZhWOt[s,:] .* Zh[:,t] ./ Zh[s,t]
      addPi += diagm(specialCaseVec)
      addPi[:,t] += specialCaseVec
      addPi[t,:] += specialCaseVec
      addPi
    end
    Pi += PlusPi
    print(t)
  end

  return(Pi)
end
