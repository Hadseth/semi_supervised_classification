############################################################
############ Sum of Similarities Classifier
############################################################
#   inputs: A - weighted adjacency matrix
#           y - all labels
#           index_train - indices of the training set
#           index_test - indices of the test set
#           alpha - diffusion parameter
#           conv_thres - threshold for convergence
#
#   ouput: node labels for the test set
############################################################

sumOfK_Classif = function(K,y,index_train,index_test)

    # Nbr of nodes
    n = length(y)
    # Nbr of classes
    nbc = maximum(y)

    # Only the training samples are labeled the others are set to a dummy class 0
    yt  = zeros(n,1)
    yt[index_train] = y[index_train]

    # Matrix n x nbc, with 1 iff node i is in class j and the training set
    Y1  = Array{Float64}(n,0)
    for c = 1:nbc
        Y1 = [Y1 (yt.==c)]
    end

    Yp = K*Y1

    # Assign to class with maximal similarity
    yp = findmax(Yp',1)[2]
    yp = (yp-1).%nbc + 1
    yp = yp'

    # Return only the label of the test set
    return(yp[index_test])
end
