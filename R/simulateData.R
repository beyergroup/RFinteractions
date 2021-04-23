simData = function(nObs,
                   nPredictors,
                   predType = c("haploid", "diploid", "continuous", "mixed"),
                   recombprob=c(0,0,0,0,0.01,0.3), # used for haploid
                   minMAF=0.3, # used for diploid and mixed
                   means=0,
                   SDs = 1, # used for continuous predictors
                   propBinary = 0.5, # used for mixed
                   outcome = c("binomial", "numeric"),
                   intType = c("none", "pure", "modifyer1", "modifyer2",
                               "redundant", "XOR", "synergistic"),
                   betaInt = c(1,1,1),
                   nMarginals,
                   betaMarginals = rep(1,nMarginals),
                   noiseSD = 1){
  # add some argument checks here.
  stopifnot("nObs must be a positive integer" = is.numeric(nObs) & nObs>0)
  # simulate predictors
  X = switch(predType,
             "haploid" = simHaploid(nObs,
                                    nPredictors,
                                    recombprob=recombprob),
             "diploid" = simDiploid(nObs,nPredictors,minMAF),
             "continuous" = simContinuous(nObs, nPredictors, means, SDs),
             "mixed" = simMixed(nObs,
                                nBinary = ceiling(propBinary*nPredictors),
                                nNumeric = floor((1-propBinary)*nPredictors),
                                means, SDs,minMAF))
  Y = simulateY(X = X,
                outcome = outcome,
                predType = predType,
                intType = intType,
                betaInt = betaInt,
                nMarginals = nMarginals,
                betaMarginals = betaMarginals,
                noiseSD = noiseSD)
  res = list(x = X, y = Y$outcome, intIDs = Y$intIDs,
             marginalIDs = Y$marginalIDs,
             intType = intType, predType = predType)
  return(res)
}

#' @title simulate haploid genotypes with tunable LD
#'
#' @description
#' A function to simulate (0,1) genotype markers with tunable linkage
#' disequilibrium
#'
#' @param nObs number of individuals
#' @param nPredictors number of predictors (i.e., markers/SNPs)
#' @param recombprob probablities of recombination; at each successive locus,
#' one of these values is sampled randomly to determine how many of the
#' individuals will "flip" their genotype; the more '0' values in this vector,
#' the more LD there will be between the markers. Higher values will result
#' in more recombination.
#' @return A dataframe with (0/1) genotypes.
#' @export
#'
#' @examples
simHaploid = function(nObs,nPredictors,
                      recombprob){
  geno = matrix(as.integer(0),nObs,nPredictors)
  xn = as.integer(sample(c(0,1),nObs,replace=T)) # initiate first marker with random genotypes
  geno[,1] = xn
  draw = ceiling(recombprob*nObs)
  for(i in 2:nPredictors){
    these = sample(1:nObs,sample(draw,1))
    xn[these] = !xn[these]
    geno[,i] = xn
  }
  return(data.frame(geno))
}


#' @title simulate diploid genotypes with tunable allele frequencies
#'
#' @description
#' A function to simulate (0,1,2) genotype markers for given allele frequencies
#'
#'
#' @param nObs number of individuals
#' @param nPredictors number of predictors (markers/SNPs)
#' @param minMAF a vector of length nPredictors with  allele frequencies. If a
#' single number, it is recycled for all predictors. If not provided (NA),
#' allele frequencies are uniformly sampled between 0 and 1.
#' @return A dataframe with (0/1/2) genotypes
#' @export
#'
#' @examples
simDiploid = function(nObs,nPredictors,minMAF){
  AFs=runif(n = nPredictors, min = minMAF, max = 1-minMAF)
  geno = sapply(AFs, function(x){
    as.integer(sample(c(0,1,2), nObs, replace = TRUE,
          prob = c((1-x)^2, 2*(1-x)*x, x^2)))
  })
  return(data.frame(geno))
}

#' @title create numeric Predictor matrix
#'
#' @description
#' A function to create a dataframe of numeric, normally distributed predictors.
#'
#'
#' @param nObs number of individuals
#' @param nPredictors number of predictors
#' @param means a vector of length nPredictors with expected mean for each
#' predictor. If a single number, it is recycled for all predictors. Defaults to 0.
#' @param SDs a vector of length nPredictors with expected standard deviation
#' for each predictor. If a single number, it is recycled for all
#' predictors. Defaults to 1.
#' @return A dataframe with numeric predictors
#' @export
#'
#' @examples
simContinuous = function(nObs, nPredictors, means=0, SDs = 1){
  if(length(means) == 1){
    means = rep(means,nPredictors)
  }
  if(length(SDs) == 1){
    SDs = rep(SDs,nPredictors)
  }
  X = sapply(seq_len(nPredictors), function(i){
    rnorm(nObs, mean = means[i], sd = SDs[i])
  })
  return(data.frame(X))
}


#' @title create numeric Predictor matrix
#'
#' @description
#' A function to create a dataframe of numeric, normally distributed predictors.
#'
#'
#' @param nObs number of individuals
#' @param nBinary number of predictors with 0/1 values.
#' @param nNumeric number of numeric/continuous predictors.
#' @param means a vector of length \code{nNumeric} with expected mean for each
#' predictor. If a single number, it is recycled for all predictors. Defaults to 0.
#' @param SDs a vector of length \code{nNumeric} with expected standard deviation
#' for each predictor. If a single number, it is recycled for all
#' predictors. Defaults to 1.
#' @param minMAF minimum allowed minor allele frequencies. If a
#' single number, it is recycled for all binary predictors If not provided (NA),
#' allele frequencies are uniformly sampled between 0 and 1.
#' @return A dataframe with \code{nBinary} numeric predictors and
#' \code{nNumeric} numeric predictors
#' @export
#'
#' @examples
simMixed = function(nObs, nBinary, nNumeric,  means=0, SDs = 1,
                    minMAF = minMAF){
  if(length(means) == 1){
    means = rep(means,nNumeric)
  }
  if(length(SDs) == 1){
    SDs = rep(SDs,nNumeric)
  }
  Xnum = sapply(seq_len(nNumeric), function(i){
    rnorm(nObs, mean = means[i], sd = SDs[i])
  })
  AFs=runif(n = nBinary, min = minMAF, max = 1-minMAF)
  Xbin = sapply(seq_len(nBinary), function(i){
    as.integer(sample(c(0,1), size = nObs, replace = T, prob = c(AFs[i],1-AFs[i])))
  })

  X = data.frame(Xnum, Xbin)
  colnames(X) = paste0("X", 1:ncol(X))
  return(X)
}



simulateY = function(X, outcome = c("binomial", "numeric"),
                     predType = c("haploid", "diploid", "continuous", "mixed"),
                     intType = c("none", "pure", "modifyer1", "modifyer2",
                                 "redundant", "XOR", "synergistic"),
                     betaInt = c(1,1,1),
                     nMarginals, betaMarginals = rep(1,nMarginals),
                     noiseSD = 1){
  # select two predictors to interact
  if(predType == "mixed"){ # for X with mixed predictor types, choose one of each
    predClasses = sapply(X,class)
    intIDs = c(sample(which(predClasses == "integer"), size = 1),
               sample(which(predClasses == "numeric"), size = 1))
  } else{
    intIDs  = sample(1:ncol(X), size = 2)
  }
  intPreds = as.matrix(X[,intIDs])
  intPreds = cbind(intPreds, intPreds[,1] * intPreds[,2])

  # select some predictors with (noise) effects
  marginalIDs = sample((1:ncol(X))[-intIDs], size = nMarginals)
  marginalPreds = as.matrix(X[,marginalIDs])
  # change the betas according to interaction type
  betaInt = switch (intType,
                    "none" = betaInt * c(1,1,0),
                    "pure" = betaInt * c(0,0,1),
                    "modifyer1" = betaInt * c(1,0,1),
                    "modifyer2" = betaInt * c(0,1,1),
                    "redundant" = betaInt * c(1,1,-1),
                    "XOR" = betaInt * c(1,1,-2),
                    "synergistic" = betaInt * c(1,1,1))
  if (outcome == "binomial"){
    zInt = mean(marginalPreds %*% cbind(betaMarginals) )
    zMarg = mean(intPreds %*% cbind(betaInt))
    # set b0 so that approximately half of cases are affected (are 1)
    b0 = -log((1-0.5)/0.5) - zInt - zMarg
    y = b0 + marginalPreds %*% cbind(betaMarginals) + intPreds %*% cbind(betaInt)
    pr = 1 / ( 1 + exp(-y))  # pass through an inv-logit function
    y = rbinom(n = nrow(X), size = 1, prob = pr) # bernoulli response variable
  } else {
    y = marginalPreds %*% cbind(betaMarginals) +
      intPreds %*% cbind(betaInt) +
      rnorm(n = nrow(X), mean = 0, sd = noiseSD)
    y = scale(y)
    y = y[,1]
  }
  return(list(outcome=y, intIDs = intIDs,
              marginalIDs = marginalIDs))
}
