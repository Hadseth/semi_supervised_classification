############################################################
############ Gaussian Kernel
############################################################


cosine_PiToK = function(Pi)
  # size
  n = size(Pi)[1]

  diagPi = diagm(1./sqrt(diag(Pi)))
  K = diagPi * Pi * diagPi
  K = (K + K')/2

  # Normalising
  Lam,U = eig(K)
  Lam,U = real(Lam),real(U)
  Lam[Lam.<0] = 0
  K = U*diagm(Lam)*U'

  return(K)
end
