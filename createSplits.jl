
createSplits = function(y,nExtFolds,nIntFolds)

  n = length(y)

  prePermSeed = randperm(n) # initial permutation
  yPerm = y[prePermSeed] # yPerm: vector y with permutation
  sortIndex = sortperm(yPerm) # index to sort vector yPerm
  permSeed = prePermSeed[sortIndex] # permSeed: index for the sorted yPerm
  extFoldsLabels = repeat(1:nExtFolds,outer=Int(ceil(n/nExtFolds)))[1:n] # Label vector for external folds
  intFoldsLabels = repeat(repeat(1:nIntFolds,inner=nExtFolds),outer=Int(ceil(n/(nExtFolds*nIntFolds))))[1:n] # Label vector for internal folds
  extFoldsLabels[permSeed] = extFoldsLabels # Label vector for external folds with initial postion
  intFoldsLabels[permSeed] = intFoldsLabels

  return(extFoldsLabels,intFoldsLabels)
end
