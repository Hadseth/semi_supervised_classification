############################################################
############ BoP Coupling
############################################################

bop_Pi = function(A,betaP;hitting=false)

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

  if hitting # Absorbing state distance (hitting paths)
        da = diag(Z)
        Za  = Z ./ (e*da')
  else # Non absorbing state distance (non hitting paths)
        Za = Z
  end

  #BoP Probability
  Pi  = Za/(e'*Za*e)

  return(Pi)
end
