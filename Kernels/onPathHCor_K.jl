############################################################
############ BoP Coupling
############################################################

onPathHCor_K = function(A, betaP)

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
      sumTot = zeros(n,n)
      for t in 1:n
        # z^{h}_{it}
        tcolZh = Zh[:,t]

        # efficient z^{h}_{it} * z^{h}_{tj}
        reshTcolZh = reshape(tcolZh,n,1)
        reshTrowZh = reshape(Zh[t,:],n,1)
        Zhittj = BLAS.gemm('N','T',reshTcolZh,reshTrowZh)

        # z^{h(-t)}_{ij}
        ZhWOt = ((Zh -  Zhittj)' ./ (1 - diag(Zhittj)))'
        ZhWOt[:,t] = 0

        # z^{h(-t)}_{\bullet i}
        rowSumZhWOt = ZhWOt'*e

        # z^{h(-t)}_{ji} * z^{h(-t)}_{ij}
        ZhWOtijji = ZhWOt.*ZhWOt'

        # z^{h(-j,-t)}_{\bullet i} * z^{h(-t)}_{ij}
        ZFiSj = ((rowSumZhWOt.*ZhWOt)' - (rowSumZhWOt.*ZhWOtijji)) ./ (1 - ZhWOtijji)
        ZFiSj[diagind(ZFiSj)] = 0

        # z^{h(i then j)}_{\bullet t}
        ZiTHENj = (ZFiSj .* tcolZh)'

        # z^{h(i then j)}_{\bullet t} + z^{h(j then i)}_{\bullet t} + \delta_{ij} z^{h(+i)}_{\bullet t}
        sumTot += ZiTHENj + ZiTHENj' + diagm(rowSumZhWOt.*tcolZh)
      end

      #### Additional terms

      #rowSum of Zh, colSum of Z, sum of Z
      rowSumZh = sum(Zh,1)[:]
      colSumZ = sum(Z,2)[:]
      sumZ = sum(Z)

      # z^{h}_{ji} * z^{h}_{ij}
      Zjiij = Zh.*Zh'

      # Z^{h(-j)}_{\bullet i} * z^{h}_{ij}
      ZbjWi = (rowSumZh.*Zh - (rowSumZh.*Zjiij)') ./ (1 .- Zjiij)
      ZbjWi[diagind(ZbjWi)] = 0

      # Eh(i,j)
      Ehij = (sumTot + ZbjWi + ZbjWi' + diagm(rowSumZh))/sumZ

      # Cov(i,j)
      ei = diag(Ehij)
      Cov = Ehij - (ei * ei')

      # Std(i)
      stdi = sqrt.(diag(Cov))

      # Cor(i and j)
      Cor = Cov ./ (stdi*stdi')

      return(Cor)
end


# onPathHCor_K = function(A,betaP)
#
#   # myMax, myMin
#   myMin = 1e-40
#   myMax = typemax(Int64)/1000
#
#   # size
#   n = size(A)[1]
#
#   # identity Matrix, e vector, centration matrix
#   I = eye(n)
#   e = ones(n)
#
#   # cost matrix, probability matrix
#   C = 1 ./ A
#   C[!isinf(C)] = C[!isinf(C)]/mean(C[!isinf(C)])
#   C[isinf(C)] = myMax
#   P = diagm(1 ./ sum(A,2)[:]) * A
#
#   # W and Z matrices
#   W = exp(-betaP * C) .* P
#   Z = inv(I - W)
#   da = diag(Z)
#
#   Zh = Z ./ (e*da')
#
#   totZh = sum(Zh)
#
#   sumTot = zeros(n,n)
#   probi = zeros(n)
#   for t in 1:n
#     Zhittj = (Zh[:,t] * Zh[t,:]')
#     ZhWOt = (Zh -  Zhittj) * diagm(1 ./ (1 - diag(Zhittj)))
#     ZhWOt[:,t] = 0
#     rowSZhWOt = sum(ZhWOt,1)[:]
#     ZFiSj = (rowSZhWOt*e' - (diagm(rowSZhWOt)*ZhWOt)') ./ (e*e' - ZhWOt.*ZhWOt')
#     ZFiSj[diagind(ZFiSj)] = 0
#     ZiTHENj = (ZFiSj .* ZhWOt) * diagm(Zh[:,t])
#     sumTot = sumTot + ZiTHENj + ZiTHENj' + diagm(rowSZhWOt.*Zh[:,t])
#     addprobi = (rowSZhWOt .* Zh[:,t])
#     addprobi[t] = sum(Zh[:,t])
#     probi = probi + addprobi
#   end
#
#   Pi = sumTot ./ totZh
#   probi = probi ./ totZh
#
#   K = Pi - (probi * probi')
#
#   varVec = sqrt(abs(diag(K)))
#   varMat = varVec * varVec'
#   K = K ./ varMat
#
#   K = (K + K')/2
#   Lam,U = eig(K)
#   Lam,U = real(Lam),real(U)
#   Lam[Lam.<0] = 0
#   K = U*diagm(Lam)*U'
#
#   return(K)
# end
