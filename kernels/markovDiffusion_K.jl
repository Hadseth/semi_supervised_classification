############################################################
############ markov diffusion Map
############################################################

# t in [2,100]

markovDiffusion_K = function(A,t)

  eps = 0.00000001

  n = size(A)[1]

  e  = ones(n)
  I  = eye(n)
  H  = (I - (e*e'/n))

  degA = sum(A,2)[:,1]
  D0   = diagm(1./ degA);
  P = D0 * A

  Z   = (inv(I - P) * (I - (P^t)) * P)/t

  K = Z * D0 * Z'
  K = H * K * H

  return(K)
end
