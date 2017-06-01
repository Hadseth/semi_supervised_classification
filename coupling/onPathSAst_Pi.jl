############################################################
############ BoP Coupling
############################################################

onPathSAst_Pi = function(A,betaP)

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
  denomZh = (1 - Zh.*Zh')
  for s in 1:n
    ZhsiMat = repmat(Zh[s,:],1,n)
    ZFiSj = ((ZhsiMat - (ZhsiMat.*Zh)') ./ denomZh)
    ZFiSj[diagind(ZFiSj)] = 0
    PlusPi = @parallel (+) for t in 1:n
      ZhjtMat = repmat(Z[:,t]',n,1)
      addPiNotSym = (ZFiSj  .* Zh .* ZhjtMat) / Z[s,t]
      addPiNotSym + addPiNotSym' + diagm(Zh[s,:] .* Z[:,t] / Z[s,t])
    end
    Pi += PlusPi
  end

  return(Pi)
end
