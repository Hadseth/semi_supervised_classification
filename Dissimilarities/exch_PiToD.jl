############################################################
######## Coupling to Dissimilarity with Exchange Matrix
############################################################

exch_PiToD = function(Pi)

  # size
  n = size(Pi)[1]

  myMax = typemax(Int64)/1000
  e = ones(n)

  PiSym = (Pi + Pi')/2
  f = sum(PiSym,2)[:,1]
  D = (diag(PiSym) ./ (f.^2)) * e' + e * (diag(PiSym) ./ (f.^2))' - 2 * (diagm(1 ./ f) * PiSym * diagm(1 ./ f) )
  D[D .> myMax/(n*(n-1))] = myMax/(n*(n-1))
  D[isinf.(D)] = myMax/(n*(n-1))

  D = D ./ mean(D)

  return(D)
end
