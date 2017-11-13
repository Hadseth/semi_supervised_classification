
nodeDistrib_D = function(A)

    # myMax, myMin
    myMin = 1e-100
    myMax = 1e+100

    # size
    n = size(A)[1]

    # identity Matrix, e vector, centration matrix
    I = eye(n)
    e = ones(n)

    ##########################

    C = 1 ./ (A + myMin)

    D = copy(C)
    for i in 1:n
        D = min(D , (D[:,i]*e' + e*D[i,:]') )
    end
    D[diagind(D)] = 0

    ##########################

    meanEdgeWeight = sum(C[A .> 0]) / sum(A .> 0)
    maxD = maximum(D)
    nbins = Int64(ceil(maxD / meanEdgeWeight))

    HistMat = zeros(n,nbins)
    for binId in 1:nbins
        countArray = (D .> (binId-1) * meanEdgeWeight) .* (D .<= binId * meanEdgeWeight)
        HistMat[:,binId] = sum(countArray,2)
    end
    HistMat = HistMat ./ (n-1)

    Ddistr = zeros(n,n)
    for i in 1:n
        xi = (HistMat[i,:]*e' + myMin)
        Ddistr[i,:] = sum(xi.*log.(2 * xi ./ (xi + HistMat')),1)
    end
    Ddistr = (Ddistr + Ddistr')/2
    Ddistr[diagind(Ddistr)] = 0

    return(Ddistr)
end
