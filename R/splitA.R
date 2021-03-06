#' @title Compute split asymmetry
#'
#' @description
#' A function to compute split asymmetry between predictors in a
#' random Forest.  For any pair of predictors (e.g. A and B), this function
#' computes the outcome variable difference (called slope) for a split on B
#' on the left side of a previous split on A, and a split on B on the right
#' side of a previous split on A. These slopes are collected over all
#' trees in the forest. The difference between left an right slopes is then
#' quantified with a Students' t-statistic (i.e. t-test).
#' Note that the p-value is not an indicator of significance for the
#' interaction, but should be treated asaan uncalibrated score, since the
#' p-values will be smaller the bigger the forest is.
#'
#' @param rf A randomForest object, created with the randomForest package.
#' @param out Either "t" or "p", indicating whether t  (out="t") or the
#' p-value (out="p") should be returned. Note that the p-value is not an
#' indicator of significance for the interaction, but should be treated as
#' an uncalibrated score.
#' @return A sparse Matrix with nPredictors*nPredictors dimensions, containing
#' either t of p-values for each predictor pair.
#' @export
#'
#' @examples
#' load("data/testRF2.RData")
#' load("data/testRF.RData")
#' res = splitA(rf)
#' res2 = splitA(rf2)
splitA = function(rf,  out = c("t", "p")){
  Mp = require(Matrix)
  RFp = require(randomForest)
  out = match.arg(out, choices = c("t", "p"))
  stopifnot("The randomForest package is required" = RFp,
            "The Matrix package is required" = Mp,
            "rf needs to be a randomForest created with randomForest package" = is(rf, "randomForest"),
            "No forest stored in rf. Use keep.forest=T when creating the forest."= is.null(rfobj$forest))

  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)
  # prepare matrix to store values in
  lSum = sparseMatrix(i = NULL, j = NULL,x=0, dims = c(nPred, nPred),
                      dimnames = list(preds, preds))
  lSqu = lSum # sum of squares of left slopes
  lCt = lSum # count of left slopes
  rSum = lSum # sum of right slopes
  rSqu = lSum # sum of squares of right slopes
  rCt = lSum # count of right slopes

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
    for(node in nodes){ # go through nodes and store all non-terminal left/right children in L/R
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
    nlSlopes = tree[tree[nl[,2],"left daughter"],"prediction"] -
      tree[tree[nl[,2],"right daughter"],"prediction"]
    nrVars = matrix(tree[nr, "split var"], ncol = 2,byrow = F)
    nrSlopes = tree[tree[nr[,2],"left daughter"],"prediction"] -
      tree[tree[nr[,2],"right daughter"],"prediction"]
    if(nrow(nlVars) < 2) print(i)
    # fill the matrixes with the slopes for this tree
    lSum[nlVars] <- lSum[nlVars] + nlSlopes
    lSqu[nlVars] <- lSqu[nlVars] + nlSlopes*nlSlopes
    lCt[nlVars] = lCt[nlVars] + 1
    rSum[nrVars] <- rSum[nrVars] + nrSlopes
    rSqu[nrVars] <- rSqu[nrVars] + nrSlopes*nrSlopes
    rCt[nrVars] = rCt[nrVars] + 1
  }


  toTest = which(lCt > 1 & rCt > 1, arr.ind = T)
  n1 = lCt[toTest]
  m1 = lSum[toTest] / n1
  v1 = (lSqu[toTest] / n1) - (m1 * m1)
  n2 = rCt[toTest]
  m2 = rSum[toTest] / n2
  v2 = (rSqu[toTest] / n2) - (m2 * m2)
  se = sqrt((((n1 - 1) * v1) + ((n2 - 1) * v2)) / (n1 + n2 - 2))
  t  = sqrt(n1 * n2 / (n1 + n2)) * ((m1 - m2) / se)

  if(out == "t"){
    res = sparseMatrix(i = toTest[,1], j = toTest[,2], x = t)
  } else{
    df <- n1 + n2 - 2
    tt_res = c(2 * pt(-abs(t), df = df))
    res = sparseMatrix(i = toTest[,1], j = toTest[,2], x = tt_res)
  }
  return(res)
}
