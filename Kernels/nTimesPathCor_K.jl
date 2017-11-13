############################################################
############ BoP Coupling
############################################################

nTimesPathCor_K = function(A,betaP)

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

    # rowSum, colSum , and sum Z
    rowSumZ = sum(Z,1)[:]
    colSumZ = sum(Z,2)[:]
    sumZ = sum(rowSumZ)

    # zbiijjb = z_{\bullet i} * z_{ij} * z_{\bullet j}
    zbiijjb = ((Z .* rowSumZ)' .* colSumZ)'

    # E(ni) = z_{\bullet i} * z_{i \bullet} / z_{\bullet \bullet}
    eni = colSumZ .* rowSumZ ./ sumZ

    # E(ni,nj)
    Eninj = (zbiijjb + zbiijjb') ./ sumZ - diagm(eni)

    # Cov(ni,nj)
    Cov = Eninj - eni*eni'

    # Std(ni)
    stdni = sqrt.(diag(Cov))

    # Cor(ni)
    Cor = Cov ./ (stdni*stdni')

    return(Cor)
end
