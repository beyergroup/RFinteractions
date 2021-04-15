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
