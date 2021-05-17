# this is a deprecated script. It is here for documentation reasons only.
######################################################################
### compile information needed for EPI() from an rfdev randomForest
### object
### rf=randomForest object (must be created with rfdev package)
### testset=logical indicating whether the testset (OOB) predictions
### or the training predictions should be used
######################################################################

#' @title simulate haploid genotypes with tunable LD
#'
#' @description
#' A function to simulate (0,1) genotype markers with tunable linkage
#' disequilibrium
#'
#' @param simData
#' @return A dataframe with (0/1) genotypes.
#'
#' @examples
EPI = function(x){
  d = abs(x$l-x$r)
  m = (abs(x$l)+abs(x$r))/2
  out = d-m
  out = out
  out = pmin(out,t(out))
  out[out<0] = 0
  return(out)
}
predsymm = function(rf){
  # create two matrices with dimensions nPredictors*nPredictors
  # to store slope sums on left and right side, respectively
  l = matrix(0,nrow=nrow(rf$importance), ncol=nrow(rf$importance))
  colnames(l) = rownames(rf$importance)
  rownames(l) = rownames(rf$importance)
  r = l
  # a matrix to count the number of collected slopes
  ct = l
  # extract infos from RF
  bvar = rf$forest$bestvar
  bvar[bvar==0]=NA
  rD = rf$forest$rightDaughter  # matrix nNodes*ntree
  rD[rD==0]=NA
  lD = rf$forest$leftDaughter # matrix nNodes*ntree
  lD[lD==0]=NA
  pred = rf$forest$nodepred # matrix nNodes*ntree

  # iterate through trees
  for(i in 1:ncol(bvar)){
    if(all(is.na(bvar[,i]))) next # if there are no splits in tree, skip this tree
    # calculate all slopes for this tree
    p = pred[lD[,i],i] - pred[rD[,i],i]
    # get the parent variable of each node
    ind = 1:nrow(bvar)
    parents = bvar[ifelse(ind %% 2 == 0, # right children will have even index
                        match(ind,lD[,i]), #match returns a vector of the positions of (first) matches of its first argument in its second.
                        match(ind,rD[,i])),i]
    # create matrix, each row a parent-child combination and the slope
    mat = cbind(nodeID = 1:nrow(bvar),parents, splitVar = bvar[,i], slope = p)
    mat = mat[apply(mat,1,function(x) all(!is.na(x))),]
    if(is.null(nrow(mat))) next
    p = mat[,4]
    ind = mat[,2:3]
    isleft = mat[,1] %% 2 != 0
    isright = !isleft
    # fill l and r with the respective slopes
    l[ind[isleft,]] = l[ind[isleft,]] + p[isleft]
    r[ind[isright,]] = r[ind[isright,]] + p[isright]
    ct[ind] = ct[ind] + 1
  }
  return(list(l=l/rf$ntree,r=r/rf$ntree,counts=ct/sum(ct)))
}
