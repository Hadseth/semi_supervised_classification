############################################################
############ BoP Coupling
############################################################

nTimesPathHCor_K = function(A,betaP)

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

  # W, Z, and Zh matrices
  W = exp.(-betaP * C) .* P
  Z = inv(I - W)
  da = diag(Z)
  Zh = Z ./ (e*da')

  # loop on t
  EninjSumTot = zeros(n,n)
  eniSumTot = zeros(n)
  for t in 1:n

      # z^{h}_{it}
      tcolZh = Zh[:,t]

      # efficient z^{h}_{it} * z_{tj}
      reshTcolZh = reshape(tcolZh,n,1)
      reshTrowZ = reshape(Z[t,:],n,1)
      Zhittj = BLAS.gemm('N','T',reshTcolZh,reshTrowZ)

      # z^{(-t)}_{ij}
      ZijWOt = Z - Zhittj

      #rowSumZijWOt = z^{(-t)}_{\bullet i}
      rowSumZijWOt = ZijWOt'*e

      # zbiijjt = z^{(-t)}_{\bullet i} * z^{(-t)}_{ij} * z_{jt}
      zbiijjt = ((ZijWOt .* rowSumZijWOt)' .* tcolZh)'

      # zbiit = z^{(-t)}_{\bullet i} * z_{it}
      zbiit = rowSumZijWOt.*tcolZh

      # symZbiijjt =  z^{(-t)}_{\bullet i} * z^{(-t)}_{ij} * z_{jt} + z^{(-t)}_{\bullet j} * z^{(-t)}_{ji} * z_{it} - \delta_{ij}
      symZbiijjt = zbiijjt + zbiijjt' - diagm(zbiit)

      EninjSumTot += symZbiijjt
      eniSumTot += zbiit
  end

  #rowSum of Z, rowSum and colSum of Zh, sum of Zh
  rowSumZ = Z'*e
  rowSumZh = Zh'*e
  sumZh = sum(rowSumZh)

  # ZbiWOjhij = z^{(-j)}_{\bullet i} * z^{h}_{ij}
  ZbiWOjhij = (rowSumZ * e' - (Z' .* rowSumZh) ) .* Zh

  # Eh(ni,nj)
  Ehninj = (EninjSumTot + ZbiWOjhij + ZbiWOjhij' + diagm(rowSumZh)) ./ sumZh

  # Eh(ni)
  ehni = (eniSumTot + rowSumZh) ./ sumZh

  # Covh(ni,nj)
  Cov = Ehninj - ehni*ehni'

  # Std(ni)
  stdhni = sqrt.(diag(Cov))

  # Corh(ni,nj)
  Cor = Cov ./ (stdhni*stdhni')

  return(Cor)
end
