# load("data/testRF2.RData")
# load("data/testRF.RData")
yates = T
selA = function(rf, yates = T, return = c("X", "p")){
  require(Matrix)
  require(randomForest)
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
  if(return == "X"){
    out = sparseMatrix(i = toTest[,1], j = toTest[,2], x = X)
  } else{
    PVAL <- pchisq(X, 1, lower.tail = FALSE)
    out = sparseMatrix(i = toTest[,1], j = toTest[,2], x = PVAL)
  }
  return(out)
}
