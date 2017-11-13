############################################################
############ Sum of Similarities Classifier
############################################################
#   inputs: A - weighted adjacency matrix
#           y - labels
#           index_train - indices of the training set
#           index_test - indices of the test set
#           alpha - diffusion parameter
#           conv_thres - threshold of convergence
#
#   ouput: node labels for the test set
############################################################

sumOfSimilarities_Classif = function(A,y,index_train,index_test,alpha;conv_thres=1e-5)

    # nbr of nodes
    n = length(y)
    # inverse out-degrees diagonal matrix
    D0 = diagm(1 ./ sum(A,2)[:])
    # transition probability matrix
    P = D0 * A
    # number of classes
    nbc = maximum(y)

    # only the training samples are labeled the others are set to a dummy class 0
    yt  = zeros(n,1)
    yt[index_train] = y[index_train]

    # matrix n x nbc, with 1 if node i is in class j and the training set
    Y1  = Array{Float64}(n,0)
    for c = 1:nbc
        Y1 = [Y1 (yt.==c)]
    end

    # Init of the degree of membership matrix n x nbc
    Yp = Y1

    # Loop
    convergence = false
    while !convergence
        Yp_prev = Yp
        # Diffision step
        Yp = alpha*P'*Yp + Y1
        # Convergence check
        if vecnorm(Yp_prev - Yp) < conv_thres
            convergence = true
        end
    end

    # Dividing results by inverse out-degrees
    Yp = D0 * Yp

    # Assign to class with maximal similarity
    yp = findmax(Yp',1)[2]
    yp = (yp-1).%nbc + 1
    yp = yp'

    # Return only the label of the test set
    return(yp[index_test])
end
