#' @title Compute paired selection frequency
#'
#' @description
#' A function to compute paired selection frequency between predictors in a
#' random Forest. This function counts, for all possible predictor pairs, how
#' often they were selected in the same tree. This count is compared to the
#' expected co-occurence, based on the individual selection frequencies, using
#' a Chi-squared statistic. Note that the p-value is not an
#' indicator of significance for the interaction, but should be treated as
#' an uncalibrated score, since the p-values will become smaller the bigger
#' the forest is.
#'
#' @param rf A randomForest object, created with the randomForest package.
#' @param yates a logical, indicating whether yates' continuity correction
#' should be applied to the Chi-squared test.
#' @param out Either "X" or "p", indicating whether Chi (out="X") or the
#' p-value (out="p") should be returned. Note that the p-value is not an
#' indicator of significance for the interaction, but should be treated as
#' an uncalibrated score.
#' @return A sparse Matrix with nPredictors*nPredictors dimensions, containing
#' either Chi-squared of p-values for each predictor pair.
#' @export
#'
#' @examples
#' load("data/testRF2.RData")
#' load("data/testRF.RData")
#' res = pairedSF(rf)
#' res2 = pairedSF(rf2)
pairedSF = function(rf, yates = T, out = c("X", "p")){
  Mp = require(Matrix)
  RFp = require(randomForest)
  out = match.arg(out, choices = c("X", "p"))
  stopifnot("The randomForest package is required" = RFp,
            "The Matrix package is required" = Mp,
            "rf needs to be a randomForest created with randomForest package" = is(rf, "randomForest"),
            "No forest stored in rf. Use keep.forest=T when creating the forest."= is.null(rfobj$forest),
            "yates needs to be a logical of length 1" = is.logical(yates) & length(yates) == 1)

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
  Y <- if (yates) { 0.5 } else { 0}
  D = abs(O-E)
  temp = cbind(Y, D)
  Y = apply(temp,1,min)
  X = rowSums(((D-Y)^2)/E)
  if(out == "X"){
    res = sparseMatrix(i = toTest[,1], j = toTest[,2], x = X)
  } else{
    PVAL <- pchisq(X, 1, lower.tail = FALSE)
    res = sparseMatrix(i = toTest[,1], j = toTest[,2], x = PVAL)
  }
  return(res)
}
