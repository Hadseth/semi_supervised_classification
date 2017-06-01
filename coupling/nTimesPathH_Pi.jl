############################################################
############ BoP Coupling
############################################################

nTimesPathH_Pi = function(A,betaP)

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
    ZWOt = Z - Zh[:,t] * Z[t,:]'
    unSym = diagm(ZWOt'*e) * ZWOt * diagm(Zh[:,t])
    vecSp = (ZWOt'*e).*(Zh[:,t])
    unSym[t,:] = unSym[t,:] + vecSp
    Sym = unSym + unSym'
    Sym[diagind(Sym)] = Sym[diagind(Sym)] - vecSp
    Sym[t,t] = Sym[t,t] + sum(Zh[:,t])
    Pi = Pi + Sym
  end

  return(Pi)
end
