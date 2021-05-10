# load("data/testRF2.RData")
# load("data/testRF.RData")

selA = function(rf){
  require(Matrix)
  require(randomForest)
  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)
  # prepare matrix to store values in
  lCt = sparseMatrix(i = NULL, j = NULL, dims = c(nPred, nPred))# sum of left slopes
  colnames(lCt) = preds
  rownames(lCt) = preds
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
    nrVars = matrix(tree[nr, "split var"], ncol = 2,byrow = F)
    # fill the matrixes with the slopes for this tree
    lCt[nlVars] = lCt[nlVars] + 1
    rCt[nrVars] = rCt[nrVars] + 1
  }
  # for each marker, how many splits were there on the left/right side in total
  nl <- rowSums(lCt)
  nr <- rowSums(rCt)

  toTest = which(lCt + rCt >= 2, arr.ind = T)

  pl = lCt/nl
  pr = rCt/nr
  # pooled Xi-squared
  p = (nl*pl+nr*pr)/(nl+nr)
  sigma = sqrt(p*(1-p)*(1/nl + 1/nr))
  # unpooled
  sigma = sqrt(pl*(1-pl)/nl + pr*(1-pr)/nr)

  z=(pl-pr)/sigma



   # prepare the counts for sel. asymmetry:
  nleft <- left[row[1],row[2],3]
  nright <- right[row[1],row[2],3]

  # binomial test:
  ntrialsL <- nTrialsLeft[row[1]]
  ntrialsR <- nTrialsRight[row[1]]

  if(ntrialsL == 0 | ntrialsR == 0){
    return(c(NA, NA, NA, NA, NA, freq_res))
  }
  chi_res <- tryCatch(unlist(prop.test(x = c(nleft,nright),
                                       n = c(ntrialsL,ntrialsR),
                                       correct = F,
                                       alternative = "two.sided")[3:4],
                              use.names = FALSE),
                       warning=function(w) c(NA, NA, NA))

}
