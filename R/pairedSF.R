# load("data/testRF2.RData")
# load("data/testRF.RData")
pairedSF = function(rf, pairMat){
  require(Matrix)
  require(randomForest)
  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)

  # prepare matrix to store values in
  ct = sparseMatrix(i = NULL, j = NULL,x=0, dims = c(nPred, nPred),
                    dimnames = list(preds, preds))
  # fill the matrix with counts
  for(i in 1:rf$ntree){
    vars = preds[unique(rf$forest$bestvar[,i])]
    ct[vars,vars] = ct[vars,vars] + 1
  }
  # only test predictor pairs that were together at least once
  toTest = which(ct > 1, arr.ind = T)
  toTest = toTest[toTest[,2]>toTest[,1],]


  # "_" stands for "not"
  n = rf$ntree
  ni = diag(ct)[toTest[,1]]
  nj = diag(ct)[toTest[,2]]
  nij = ct[toTest]
  ni_j = ni - nij
  n_ij = nj - nij
  n_i_j = n - ni - nj + nij
  n_i = n-ni
  n_j = n-nj
  O = cbind(nij, ni_j, n_ij, n_i_j) # observed values
  E = cbind(ni*nj,ni*n_j,n_i*nj, n_i*n_j)/n # expected values
  X = rowSums(((O-E)^2)/E)
  return(X)
}
