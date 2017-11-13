############################################################
############ boPLabelProp_Classif
############################################################

BoPLabelProp_Classif = function(A,y,indexTrain,indexTest,betaP;hitting = false,rev=false)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(A)[1]

  # number of class
  nbc = maximum(y)

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)

  # supply and demand
  sIn = zeros(n)
  sOut = zeros(n)
  if rev
    sIn[indexTest] = 1
    sOut[indexTrain] = 1
  else
    sIn[indexTrain] = 1
    sOut[indexTest] = 1
  end

  # cost matrix, probability matrix
  C = 1 ./ A
  C[.!isinf.(C)] = C[.!isinf.(C)]/mean(C[.!isinf.(C)])
  C[isinf.(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:]) * A

  # W and Z matrices
  W = exp.(-betaP * C) .* P
  Z = inv(I - W)

  if hitting # Absorbing state distance (hitting paths)
    da = diag(Z)
    Za  = Z ./ (e*da')
  else # Non absorbing state distance (non hitting paths)
    Za = Z
  end

  #BoP Probability
  Pi  = Za/(e'*Za*e)

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
