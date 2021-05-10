# load("data/testRF2.RData")
# load("data/testRF.RData")
pairedSF = function(rf, pairMat){
  require(Matrix)
  require(randomForest)
  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)
  
  
  # prepare matrix to store values in
  ct =  sparseMatrix(i = NULL, j = NULL,x=0, dims = c(nPred, nPred),
                     dimnames = list(preds, preds))
 
  for(i in 1:rf$ntree){
    vars = preds[unique(rf$forest$bestvar[,i])]
    ct[vars,vars] = ct[vars,vars]+ 1
    # ct[-vars,vars,2] = ct[-markers,markers,2] + 1
  }
  
  # increase the counter for the marker pairs that occured together
  # counts[markers,markers,1] = counts[markers,markers,1] + 1
  # # increase the counter for the cases where only one of the markers was used (not used in row)
  # counts[-markers,markers,2] = counts[-markers,markers,2] + 1
  
  
  # diag(ct) = 0
  toTest = which(ct > 1, arr.ind = T)
  i = 1
  j = 3
  nij = ct[i,j]
  ni = ct[i,i]-nij
  nj = ct[j,j]-nij
  nr = rf$ntree - ni - nj - nij
  mat <- matrix(c(nij, ni, nj, nr), nrow = 2)
  chisq.test(mat)$statistic
}