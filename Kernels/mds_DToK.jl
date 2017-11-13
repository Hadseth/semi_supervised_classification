############################################################
############ MDS Kernel
############################################################

mds_DToK = function(D)

  # size
  n = size(D)[1]

  # identity Matrix, e vector, centration matrix
  I = eye(n)
  e = ones(n)
  H = I - e*e' / n

  # Compute the MDS Kernel
  K = -0.5 * H * (D.^2) * H

  # Normalising
  Lam,U = eig(K)
  Lam,U = real(Lam),real(U)
  Lam[Lam.<0] = 0
  K = U*diagm(Lam)*U'

  return(K)
end
