############################################################
############ BoP Coupling
############################################################

freeEnergy_D = function(A,betaP)

  # myMax, myMin
  myMin = 1e-40
  myMax = typemax(Int64)/1000

  # size
  n = size(A)[1]

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)
  H = I - e*e' / n

  # cost matrix, probability matrix
  C = 1 ./ A
  C[!isinf(C)] = C[!isinf(C)]/mean(C[!isinf(C)])
  C[isinf(C)] = myMax
  P = diagm(1 ./ sum(A,2)[:,1]) * A

  # W and Z matrices
  W = exp(-betaP * C) .* P
  Z = inv(I - W)

  da = diag(Z)
  Za  = Z ./ (e*da')

  # Dissimilarity Construction
  Zpart = Za
  Zpart[Zpart .< myMin] = myMin
  D = -log(Zpart)/betaP
  D = 0.5 * (D + D')
  D = D - diagm(diag(D))

  D = D ./ mean(D)

  return(D)
end
