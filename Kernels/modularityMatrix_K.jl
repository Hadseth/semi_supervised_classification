############################################################
############ modularity Matrix Kernel
############################################################

modularityMatrix_K = function(A)
  # size
  n = size(A)[1]

  # graph volume
  vol = sum(A)

  # in and out degree
  din  = sum(A,1)
  dout = sum(A,2)

  D  = (dout * din)/vol

  K  = (A - D)

  return(K)
end
