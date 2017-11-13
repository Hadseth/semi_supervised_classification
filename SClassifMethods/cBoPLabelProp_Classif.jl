############################################################
############ boPLabelProp_Classif
############################################################

cBoPLabelProp_Classif = function(A,y,indexTrain,indexTest,betaP;hitting = false,sInW=repeat([1/(size(A)[1])],inner=size(A)[1]),sOutW=repeat([1/(size(A)[1])],inner=size(A)[1]),rev=false,persist=1,convThres=1e-10)

    # myMax, myMin
    myMin   = 1e-40
    myMax = typemax(Int64)/1000

    # size
    n = size(A)[1]

    # number of class
    nbc = maximum(y)

    # identity Matrix, e vector, centration matrix
    I = eye(n)
    e = ones(n)

    # cost matrix, probability matrix
    C = 1 ./ A
    C[.!isinf.(C)] = C[.!isinf.(C)]/mean(C[.!isinf.(C)])
    C[isinf.(C)] = myMax
    P = diagm(1 ./ sum(A,2)[:,1]) * A

    sIn = zeros(n)
    sOut = zeros(n)
    if rev
      sIn[indexTest] = sInW[indexTest]
      sOut[indexTrain] = sOutW[indexTrain]
    else
      sIn[indexTrain] = sInW[indexTrain]
      sOut[indexTest] = sOutW[indexTest]
    end

    sIn = sIn / sum(sIn)
    sOut = sOut / sum(sOut)

    if hitting
      W = exp.(-betaP * C) .* P
      Z = inv(I - W)
      Z = Z * diagm(1 ./ diag(Z))
      outVal = sOut
      numerator = e
    else
      Q = I - P'
      pinvQ = pinv(Q)
      q = (I - pinvQ * Q)
      q = q[:,1]

      nRef0 = pinvQ * (sIn - P' * sOut)
      gammaP = maximum((sOut - nRef0) ./ q) + persist
      nRef = nRef0 + gammaP * q
      alpha = sOut ./ nRef
      Pmod = (I - diagm(alpha)) * P

      # matrices W and Z, vectors divisor and out_val
      W = exp.(-betaP * C) .* Pmod
      Z = inv(I - W)
      outVal = alpha
      numerator = nRef
    end

    # Iterations
    mIn = e
    mOut = e
    convergence = false
    FE = myMax
    while !convergence
      FEPrev = FE
      mIn = e ./ (Z * (mOut .* outVal))
      mOut = numerator ./ (Z' * (mIn .* sIn))
      FE = sum(-log.(mIn) .* sIn) + sum(-log.(mOut) .* sOut)
      maxDiff = abs(FEPrev - FE)
      if maxDiff < convThres
        convergence = true
      end
    end
    Pi = diagm(mIn .* sIn) * Z * diagm(mOut .* outVal)

    if rev
      sumByClass = Array{Float64}(0,sum((sIn .> 0)))
      for class in 1:nbc
        sumByClass = [sumByClass;sum(Pi[indexTest,((sOut .> 0) & (y .== class))[:]],2)']
      end
    else
      sumByClass = Array{Float64}(0,sum((sOut .> 0)))
      for class in 1:nbc
        sumByClass = [sumByClass;sum(Pi[((sIn .> 0) & (y .== class))[:],indexTest],1)]
      end
    end
    yp = findmax(sumByClass,1)[2]
    labelPrediction = ((yp-1).%nbc + 1)[:]

  return(labelPrediction)
end
