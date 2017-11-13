############################################################
############ BoP Coupling
############################################################

onPathCov_K = function(A,betaP)

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
  C[.!isinf.(C)] = C[.!isinf.(C)]/mean(C[.!isinf.(C)])
  C[isinf.(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:]) * A

  # W and Z matrices
  W = exp.(-betaP * C) .* P
  Z = inv(I - W)
  da = diag(Z)

  Zh = Z ./ (e*da')

  ######## Covariance

  #rowSum of Zh, colSum of Z, sum of Z
  rowSumZh = sum(Zh,1)[:]
  colSumZ = sum(Z,2)[:]
  sumZ = sum(Z)

  # z^{h}_{ji} * z^{h}_{ij}
  Zjiij = Zh.*Zh'

  # Z^{h(+i)}_{\bullet j} = Z^{h(-j)}_{\bullet i} * z^{h}_{ij} + delta_{ij} * z^{h}_{\bullet i}
  ZbjWi = (rowSumZh.*Zh - (rowSumZh.*Zjiij)') ./ (1 .- Zjiij)
  ZbjWi[diagind(ZbjWi)] = rowSumZh

  # E(i then j)
  EijNotSym = (ZbjWi' .* colSumZ)' / sumZ

  # E(i)
  ei = diag(EijNotSym)

  # E(i and j)
  Eij = EijNotSym + EijNotSym' - diagm(ei)

  # Cov(i and j)
  Cov =  Eij - ei*ei'

  return(Cov)
end
