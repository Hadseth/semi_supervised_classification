############################################################
############ MDS Kernel
############################################################

invMDS_KToD = function(K)

  # size
  n = size(K)[1]

  # e vector
  e = ones(n)

  # diag K
  diagK = diag(K)

  # D
  D = diagK*e' + e*diagK' - (2 .* K)

  # securities
  D[D .< 0] = 0
  D = sqrt.(D)
  D[diagind(D)] = 0
  D = (D + D')/2

  return(D)

end
