############################################################
############ regularized Commute Time Kernel
############################################################

# a  = [0.01,0.99]

regularizedCommuteTime_K = function(A,a)
  eps = 0.01

  if a < eps
    a = eps
  elseif a > (1-eps)
    a = (1-eps)
  end

  # size
  n = size(A)[1]

  I = eye(n)
  D  = diagm(sum(A,2)[:])

  K  = inv(D - (a*A))

  return(K)
end
