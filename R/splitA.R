# load("data/testRF2.RData")
# load("data/testRF.RData")


splitA = function(rf){
  require(Matrix)
  require(randomForest)
  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)
  # prepare matrix to store values in
  lSum = sparseMatrix(i = NULL, j = NULL, dims = c(nPred, nPred))# sum of left slopes
  colnames(lSum) = preds
  rownames(lSum) = preds
  lSqu = lSum # sum of squares of left slopes
  lCt = lSum # count of left slopes
  rSum = lSum # sum of right slopes
  rSqu = lSum # sum of squares of right slopes
  rCt = lSum # count of right slopes

  #organize the forest-matrices into an array
  # forestArray <- array(0,dim = c(dim(rf$forest$nodestatus),5))
  # forestArray[,,1] <- rf$forest$leftDaughter # rows: node, cols: tree
  # forestArray[,,2] <- rf$forest$rightDaughter # rows: node, cols: tree
  # forestArray[,,3] <- rf$forest$nodestatus # whether it is a terminal node or not (-3 means not terminal, -1 terminal, and 0 not used)
  # forestArray[,,4] <- rf$forest$bestvar # index of the alleles that were used for the split at each node (and in each tree)
  # forestArray[,,5] <- rf$forest$oobpred

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
  df <- n1 + n2 - 2
  tt_res = c(2 * pt(-abs(t), df = df))
  out = sparseMatrix(i = toTest[,1], j = toTest[,2], x = tt_res) # use t instead of tt_res?
  return(out)
}
