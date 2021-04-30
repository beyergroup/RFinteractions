# load("data/testRF2.RData")
# load("data/testRF.RData")

# function to extract Jake's score from a rf. Improved to also include
# scores for predictors that are not directly consecutie in the trees.
slopeScore = function(rf){
  require(Matrix)
  require(randomForest)
  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)
  l = sparseMatrix(i = NULL, j = NULL, dims = c(nPred, nPred))
  ## the (i,j) pairs can be repeated, in which case the x's are summed
  # l2=   matrix(0,nrow=nrow(rf$importance), ncol=nrow(rf$importance)) # empty matrix is much, much bigger
  colnames(l) = preds
  rownames(l) = preds
  r = l
  # matrices to count the number of collected slopes
  ctL = l
  ctR = l
  for(i in 1:rf$ntree){
    # extract tree information from forest
    if (rf$type == "regression") {
      tree <- cbind(rf$forest$leftDaughter[,i],
                    rf$forest$rightDaughter[,i],
                    rf$forest$bestvar[,i],
                    rf$forest$xbestsplit[,i],
                    rf$forest$nodestatus[,i] == -3,
                    rf$forest$nodepred[,i])[1:rf$forest$ndbigtree[i], ]
    }else {
      tree <- cbind(rf$forest$treemap[, , i],
                    rf$forest$bestvar[, i],
                    rf$forest$xbestsplit[, i],
                    rf$forest$nodestatus[, i] == 1,
                    rf$forest$nodepred[, i])[1:rf$forest$ndbigtree[i], ]
    }
    dimnames(tree) <- list(1:nrow(tree),
                           c("left daughter","right daughter",
                             "split var", "split point",
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
        L[[node]] <- c(as.integer(dirChildren[1]), L[[dirChildren[1]]], R[[dirChildren[1]]])
      }
      if(logic[2]){
        R[[node]] <- c(as.integer(dirChildren[2]), L[[dirChildren[2]]], R[[dirChildren[2]]])
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
    #build a matrix: col 1&2: nodes, col3&4: respective predictors
    nlVars = matrix(tree[nl, "split var"], ncol = 2,byrow = F)
    nlSlopes = tree[tree[nl[,2],"left daughter"],"prediction"] -
      tree[tree[nl[,2],"right daughter"],"prediction"]
    nrVars = matrix(tree[nr, "split var"], ncol = 2,byrow = F)
    nrSlopes = tree[tree[nr[,2],"left daughter"],"prediction"] -
      tree[tree[nr[,2],"right daughter"],"prediction"]
    l[nlVars] <- l[nlVars] + nlSlopes
    ctL[nlVars] = ctL[nlVars] + 1
    r[nrVars] <- r[nrVars] + nrSlopes
    ctR[nlVars] = ctR[nlVars] + 1
  }
  d = abs(l-r)
  m = (abs(l)+abs(r))/2
  out = d-m
  # out = pmin(out,t(out))
  out[out<0] = 0
  return(out)
}
