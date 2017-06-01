############################################################
############ BoP Coupling
############################################################

onPath_Pi = function(A,betaP)

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

  ZFiSj = (Zh'*e*e' - (diagm(Zh'*e)*Zh)') ./ (e*e' - Zh.*Zh')
  ZFiSj[diagind(ZFiSj)] = 0
  PiNotSym = (ZFiSj .* Zh) * diagm(Z*e)

  #Probability to have i and j on non-hitting paths
  Pi  = PiNotSym + PiNotSym'
  Pi = Pi + diagm((Zh'*e) .* (Z*e))

  Pi = Pi ./ sum(Z)

  return(Pi)
end
