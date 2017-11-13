############################################################
############ regularized Laplacian Kernel
############################################################

# a  = [0.001,0.999]

regularizedLaplacian_K = function(A,a)
  eps = 0.001

  if a < eps
    a = eps
  elseif a > (1-eps)
    a = (1-eps)
  end

  # size
  n = size(A)[1]

  I = eye(n)
  D  = diagm(sum(A,2)[:])
  L = D - A

  K  = inv(I + (a*L))

  return(K)
end
