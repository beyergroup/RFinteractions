#' @title Compute slope asymmetry score
#'
#' @description
#' A function to compute aggregated slope asymmetry scores between predictors in a
#' random Forest.  For any pair of predictors (e.g. A and B), this function
#' computes the outcome variable difference (called slope) for a split on B
#' on the left side of a previous split on A, and a split on B on the right
#' side of a previous split on A. These slopes are aggregated (summed) over all
#' trees in the forest. The difference between left an right slopes is then
#' returned, after correcting for marginal effects of each predictor.
#' This function represents an improved implementation of the
#' approach presented in the PhD Thesis of Jake Michaelson. The major
#' modification is also including predictor pairs that are not directly
#' consecutive in the tree, but also predictor pairs where there were splits
#' on other predictors in between.
#' @param rf A randomForest object, created with the randomForest package.
#' @return A sparse Matrix with nPredictors*nPredictors dimensions, containing
#' slope asymmetry scores for predictor pairs. Rows represent parent split (A),
#' columns represent the child split (B). Therefore, each predictor pair is
#' represented twice, once with A as the parent, and once with B as the parent.
#' @export
#'
#' @examples
#' load("data/testRF2.RData")
#' load("data/testRF.RData")
#' res = selA(rf)
#' res2 = selA(rf2)
slopeScore = function(rf){
  Mp = require(Matrix)
  RFp = require(randomForest)
  stopifnot("The randomForest package is required" = RFp,
            "The Matrix package is required" = Mp,
            "rf needs to be a randomForest created with randomForest package" = is(rf, "randomForest"),
            "No forest stored in rf. Use keep.forest=T when creating the forest."= is.null(rfobj$forest))

  nPred = nrow(rf$importance)
  preds = rownames(rf$importance)
  l =  sparseMatrix(i = NULL, j = NULL,x=0, dims = c(nPred, nPred),
                    dimnames = list(preds, preds))
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
