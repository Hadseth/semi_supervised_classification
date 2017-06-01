############################################################
############ BoP Coupling
############################################################

onPathH_Pi = function(A,betaP)

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

  sumTot = zeros(n,n)
  for t in 1:n
    Zhittj = (Zh[:,t] * Zh[t,:]')
    ZhWOt = (Zh -  Zhittj) * diagm(1 ./ (1 - diag(Zhittj)))
    ZhWOt[:,t] = 0
    rowSZhWOt = sum(ZhWOt,1)[:]
    ZFiSj = (rowSZhWOt*e' - (diagm(rowSZhWOt)*ZhWOt)') ./ (e*e' - ZhWOt.*ZhWOt')
    ZFiSj[diagind(ZFiSj)] = 0
    ZiTHENj = (ZFiSj .* ZhWOt) * diagm(Zh[:,t])
    sumTot = sumTot + ZiTHENj + ZiTHENj' + diagm(rowSZhWOt.*Zh[:,t])
  end

  Pi = sumTot ./ sum(Zh)
  return(Pi)
end
