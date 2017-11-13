############################################################
############ Gaussian Kernel
############################################################


gauss_DToK = function(D)
  # size
  n = size(D)[1]

  # Compute the Gaussian Kernel
  sig = sum(D)/(n*(n-1))
  K = exp.(- D ./sig)

  # Normalising
  Lam,U = eig(K)
  Lam,U = real(Lam),real(U)
  Lam[Lam.<0] = 0
  K = U*diagm(Lam)*U'

  return(K)
end
