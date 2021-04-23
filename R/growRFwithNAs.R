#' @title Grow a Random Forest for data with missing values
#' 
#' @description
#' A wrapper for RandomForestExtended::randomForest that takes care of the handling
#' of missing values. 
#' When there are missing values in the outcome variable (y),
#' they are removed from x and y. Missing values in the predictors (x) are 
#' imputed as follows:
#' For each of \code{nImpute} iterations, missing values in x were replaced 
#' with values that were randomly sampled (with replacement) from the same column
#' of x (i.e. the same predictor). A small sub-Random Forest is grown on this imputed data.
#' Finally, all \code{nImput} sub-Random Forests are combined into one big 
#' Random Forest.
#' 
#' This approach is essentially equivalent to using random split assignments for 
#' missing predictors. There higher nImpute, the lower the likelyhood that 
#' a single imputation iteraction has a large impact on the final RF.
#' 
#' @param x A data frame or a matrix of predictors, or a formula describing 
#' the model to be fitted (for the print method, an randomForest object).
#' If x is a matrix it can contain either double or raw values. 
#' The use of raw values needs much less memory.
#' @param y A response vector. If a factor, classification is assumed, 
#' otherwise regression is assumed. If omitted, randomForest will 
#' run in unsupervised mode. 
#' @param ntree number of trees in the final Random Forest
#' @param nImpute A positive integer indicating how often data should be
#' imputed. It is recommended that ntree is a multiple of nImpute, otherwise
#' the final forest size might differ from ntree, with a warning.
#' @param ... passed on to \code{randomForest}.
#' @return An object of class randomForest.
#' @export
#'
#' @examples
growRFwithNAs = function(x, y, nImpute=10, ntree=500, ...){
  RFE = require(RandomForestExtended)
  stopifnot("RFepistasis depends on the package RandomForestExtended. 
           It can be accessed at http://cellnet-sb.cecad.uni-koeln.de/resources/RandomForestExtended." = RFE,
          "nImpute has to be a positive integer." = nImpute<1 | !is.numeric(nImpute),)
  if(ntree %% nImpute != 0)
    warning("ntree is not a multiple of nImpute:\n final RF size may differ from ntree")
  # all other arguments are checked within the randomForest function
  missingY = is.na(y)
  y = y[!missingY]
  rfList = lapply(1:nImpute, FUN = function(x){
    xImputed = sapply(1:ncol(x),FUN=function(i){ #per predictor
      missingX = is.na(i)
      i[missingX] = sample(i[!missingX],size = sum(missingX), replace = T)
      return(i) #return the extended genotype
    })
    xImputed = xImputed[!missingY,]
    RandomForestExtended::randomForest(x=xImputed,
                                       y=y,
                                       ntree=ceiling(ntree/nImpute),
                                       ...)
  })
  rf = combineForestsExtended(rfList)
  return(rf)
}

#' @title Combine RFs
#' 
#' @description
#' Function that combines several randomForest objects 
#' into one big Random Forest. Workaround for a bug in the
#' function  \code{combine} in \code{RandomForestExtended}.
#'
#' @param ... Two or more objects of class randomForest, to be combined into one.
#' 
#' @return An object of class randomForest.
#' @export
#'
#' @examples
combineForestsExtended = function (...) {
  pad0 <- function(x, len) c(x, rep(0, len - length(x)))
  padm0 <- function(x, len) rbind(x, matrix(0, nrow = len - 
                                              nrow(x), ncol = ncol(x)))
  rflist <- list(...)
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) 
    rflist <- rflist[[1]]
  areForest <- sapply(rflist, function(x) inherits(x, "randomForest"))
  if (any(!areForest)) 
    stop("Argument must be a list of randomForest objects")
  rf <- rflist[[1]]
  classRF <- rf$type == "classification"
  trees <- sapply(rflist, function(x) x$ntree)
  ntree <- sum(trees)
  rf$ntree <- ntree
  nforest <- length(rflist)
  haveTest <- !any(sapply(rflist, function(x) is.null(x$test)))
  vlist <- lapply(rflist, function(x) rownames(RandomForestExtended::importance(x)))
  numvars <- sapply(vlist, length)
  if (!all(numvars[1] == numvars[-1])) 
    stop("Unequal number of predictor variables in the randomForest objects.")
  for (i in seq_along(vlist)) {
    if (!all(vlist[[i]] == vlist[[1]])) 
      stop("Predictor variables are different in the randomForest objects.")
  }
  haveForest <- sapply(rflist, function(x) !is.null(x$forest))
  if (all(haveForest)) {
    nrnodes <- max(sapply(rflist, function(x) x$forest$nrnodes))
    rf$forest$nrnodes <- nrnodes
    rf$forest$ndbigtree <- unlist(sapply(rflist, function(x) x$forest$ndbigtree))
    rf$forest$nodestatus <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodestatus, nrnodes)))
    rf$forest$bestvar <- do.call("cbind", lapply(rflist, 
                                                 function(x) padm0(x$forest$bestvar, nrnodes)))
    rf$forest$xbestsplit <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$xbestsplit, nrnodes)))
    rf$forest$nodepred <- do.call("cbind", lapply(rflist, 
                                                  function(x) padm0(x$forest$nodepred, nrnodes)))
    tree.dim <- dim(rf$forest$treemap)
    if (classRF) {
      rf$forest$treemap <- array(unlist(lapply(rflist, 
                                               function(x) apply(x$forest$treemap, 2:3, pad0, 
                                                                 nrnodes))), c(nrnodes, 2, ntree))
    }
    else {
      rf$forest$leftDaughter <- do.call("cbind", lapply(rflist, 
                                                        function(x) padm0(x$forest$leftDaughter, nrnodes)))
      rf$forest$rightDaughter <- do.call("cbind", lapply(rflist, 
                                                         function(x) padm0(x$forest$rightDaughter, nrnodes)))
      rf$forest$nodepops <- do.call("cbind", lapply(rflist, 
                                                    function(x) padm0(x$forest$nodepops, nrnodes)))
      rf$forest$oobpred <- do.call("cbind", lapply(rflist, 
                                                   function(x) padm0(x$forest$oobpred, nrnodes)))
      rf$forest$oobmse <- do.call("cbind", lapply(rflist, function(x) padm0(x$forest$oobmse, 
                                                                            nrnodes)))
    }
    rf$forest$ntree <- ntree
    if (classRF) 
      rf$forest$cutoff <- rflist[[1]]$forest$cutoff
  }
  else {
    rf$forest <- NULL
  }
  if (classRF) {
    rf$votes <- 0
    rf$oob.times <- 0
    areVotes <- all(sapply(rflist, function(x) any(x$votes > 
                                                     1, na.rf = TRUE)))
    if (areVotes) {
      for (i in 1:nforest) {
        rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
        rf$votes <- rf$votes + ifelse(is.na(rflist[[i]]$votes), 
                                      0, rflist[[i]]$votes)
      }
    }
    else {
      for (i in 1:nforest) {
        rf$oob.times <- rf$oob.times + rflist[[i]]$oob.times
        rf$votes <- rf$votes + ifelse(is.na(rflist[[i]]$votes), 
                                      0, rflist[[i]]$votes) * rflist[[i]]$oob.times
      }
      rf$votes <- rf$votes/rf$oob.times
    }
    rf$predicted <- factor(colnames(rf$votes)[max.col(rf$votes)], 
                           levels = levels(rf$predicted))
    if (haveTest) {
      rf$test$votes <- 0
      if (any(rf$test$votes > 1)) {
        for (i in 1:nforest) rf$test$votes <- rf$test$votes + 
            rflist[[i]]$test$votes
      }
      else {
        for (i in 1:nforest) rf$test$votes <- rf$test$votes + 
            rflist[[i]]$test$votes * rflist[[i]]$ntree
      }
      rf$test$predicted <- factor(colnames(rf$test$votes)[max.col(rf$test$votes)], 
                                  levels = levels(rf$test$predicted))
    }
  }
  else {
    rf$predicted <- 0
    for (i in 1:nforest) rf$predicted <- rf$predicted + rflist[[i]]$predicted * 
        rflist[[i]]$ntree
    rf$predicted <- rf$predicted/ntree
    if (haveTest) {
      rf$test$predicted <- 0
      for (i in 1:nforest) rf$test$predicted <- rf$test$predicted + 
          rflist[[i]]$test$predicted * rflist[[i]]$ntree
      rf$test$predicted <- rf$test$predicted/ntree
    }
  }
  have.imp <- !any(sapply(rflist, function(x) is.null(x$importance)))
  if (have.imp) {
    rf$importance <- rf$importanceSD <- 0
    for (i in 1:nforest) {
      rf$importance <- rf$importance + rflist[[i]]$importance * 
        rflist[[i]]$ntree
      rf$importanceSD <- rf$importanceSD + rflist[[i]]$importanceSD^2 * 
        rflist[[i]]$ntree
    }
    rf$importance <- rf$importance/ntree
    rf$importanceSD <- sqrt(rf$importanceSD/ntree)
    haveCaseImp <- !any(sapply(rflist, function(x) is.null(x$localImportance)))
    if (haveCaseImp) {
      rf$localImportance <- 0
      for (i in 1:nforest) {
        rf$localImportance <- rf$localImportance + rflist[[i]]$localImportance * 
          rflist[[i]]$ntree
      }
      rf$localImportance <- rf$localImportance/ntree
    }
  }
  have.prox <- !any(sapply(rflist, function(x) is.null(x$proximity)))
  if (have.prox) {
    rf$proximity <- 0
    for (i in 1:nforest) rf$proximity <- rf$proximity + rflist[[i]]$proximity * 
        rflist[[i]]$ntree
    rf$proximity <- rf$proximity/ntree
  }
  if (classRF) {
    rf$confusion <- NULL
    rf$err.rate <- NULL
    if (haveTest) {
      rf$test$confusion <- NULL
      rf$err.rate <- NULL
    }
  }
  else {
    rf$mse <- rf$rsq <- NULL
    if (haveTest) 
      rf$test$mse <- rf$test$rsq <- NULL
  }
  keep.inbag <- !is.null(rf$inbag)
  if (keep.inbag) {
    for (i in 2:length(rflist)) {
      rf$inbag <- cbind(rf$inbag, rflist[[i]]$inbag)
    }
  }
  else {
    rf$inbag <- NULL
  }
  rf
}

