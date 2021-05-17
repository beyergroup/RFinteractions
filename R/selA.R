#' @title Compute selection asymmetry
#'
#' @description
#' A function to selection asymmetry between predictors in a
#' random Forest. This function counts, for all possible predictor pairs, how
#' often on predictor was used on the left or the right side of a previous
#' split on the second predictor. The imbalance in these counts is quantified
#' with a Chi-squared statistic. Note that the p-value is not an
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
#' res = selA(rf)
#' res2 = selA(rf2)
selA = function(rf, yates = T, out = c("X", "p")){
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
  lCt =  sparseMatrix(i = NULL, j = NULL,x=0, dims = c(nPred, nPred),
                      dimnames = list(preds, preds))
  rCt = lCt # count of right slopes

  for(i in 1:rf$ntree){
    # extract tree information from forest
    if (rf$type == "regression") {
      tree <- cbind(rf$forest$leftDaughter[,i],
                    rf$forest$rightDaughter[,i],
                    rf$forest$bestvar[,i],
                    rf$forest$nodestatus[,i] == -3,
                    rf$forest$nodepred[,i])[1:rf$forest$ndbigtree[i], ]
    }else {
      tree <- cbind(rf$forest$treemap[, , i],
                    rf$forest$bestvar[, i],
                    rf$forest$nodestatus[, i] == 1,
                    rf$forest$nodepred[, i])[1:rf$forest$ndbigtree[i], ]
    }
    dimnames(tree) <- list(1:nrow(tree),
                           c("left daughter","right daughter",
                             "split var",
                             "status", "prediction"))
    # get only nodes which have children
    subTree <- tree[tree[,"status"] == 1,, drop = F]
    # create lists L and R with recursive left and right children
    L <- vector("list",nrow(subTree))
    names(L) <- rownames(subTree)
    R <- L
    nodes <- rev(rownames(subTree))
    # go through nodes and store all non-terminal left/right children in L/R
    for(node in nodes){
      dirChildren <- as.character.default(subTree[node,1:2])
      logic <- (dirChildren %in% nodes) # does it have non-terminal left or right children?
      if(logic[1]){
        L[[node]] <- c(as.integer(dirChildren[1]),
                       L[[dirChildren[1]]],
                       R[[dirChildren[1]]])
      }
      if(logic[2]){
        R[[node]] <- c(as.integer(dirChildren[2]),
                       L[[dirChildren[2]]],
                       R[[dirChildren[2]]])
      }
    }

    ## transform above lists into matrices with valid marker combinations,
    ##  first column mother, second column (direct or indirect) daughter node
    nl <- c(unlist(sapply(names(L), FUN = function(x){
      rep(as.numeric(x), length(L[[x]]))
    },USE.NAMES = F)), unlist(L, use.names = FALSE))
    if(is.null(nl)){
      nl = matrix(data = NA, nrow = 0, ncol = 2)
    } else{
      dim(nl) <- c(length(nl)/2, 2)
    }
    nr <- c(unlist(sapply(names(R), FUN = function(x){
      rep(as.numeric(x), length(R[[x]]))
    },USE.NAMES = F)), unlist(R, use.names = FALSE))
    if(is.null(nr)){
      nr = matrix(data = NA, nrow = 0, ncol = 2)
    }else{
      dim(nr) <- c(length(nr)/2, 2)
    }
    # get indices and slopes for predictor pairs
    nlVars = matrix(tree[nl, "split var"], ncol = 2,byrow = F)
    nrVars = matrix(tree[nr, "split var"], ncol = 2,byrow = F)
    # fill the matrixes with the slopes for this tree
    lCt[nlVars] = lCt[nlVars] + 1
    rCt[nrVars] = rCt[nrVars] + 1
  }

  # which marker pairs to Test
  toTest = which(lCt + rCt >= 2, arr.ind = T)
  # we mustn't do the following:
  # toTest = toTest[toTest[,2]>toTest[,1],]

  nTrialsL = rowSums(lCt)[toTest[,1]]
  nTrialsR = rowSums(rCt)[toTest[,1]]
  nl = lCt[toTest]
  nr = rCt[toTest]
  Y <- if (yates) { 0.5 } else { 0}
  delta = abs(nl/nTrialsL - nr/nTrialsR)
  Y <- pmin(Y, delta / (1/nTrialsL + 1/nTrialsR))
  p = (nl+nr)/(nTrialsL + nTrialsR)
  n_l = nTrialsL - nl
  n_r = nTrialsR - nr
  O = cbind(nl, nr, n_l, n_r)
  E = cbind(p*nTrialsL, p*nTrialsR, (1-p) * nTrialsL, (1-p) * nTrialsR)
  X <- rowSums((abs(O - E) - Y)^2/E)
  if(out == "X"){
    resetGeneric() = sparseMatrix(i = toTest[,1], j = toTest[,2], x = X)
  } else{
    PVAL <- pchisq(X, 1, lower.tail = FALSE)
    resetClass() = sparseMatrix(i = toTest[,1], j = toTest[,2], x = PVAL)
  }
  return(res)
}
