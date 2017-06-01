############################################################
############ Surprisal Dissimilarity from coupling
############################################################

sur_PiToD = function(Pi)

  # size
  n = size(Pi)[1]

  myMax = typemax(Int64)/1000

  D = -(log(Pi) + log(Pi'))/2
  D = D- diagm(diag(D))
  D[D .> myMax/(n*(n-1))] = myMax/(n*(n-1))
  D[isinf(D)] = myMax/(n*(n-1))

  D = D ./ mean(D)
  return(D)
end
