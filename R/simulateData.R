#' @title simulate data for the benchmark of interaction detection methods
#'
#' @description
#' A function to simulate different types of variables as predictors and their
#' influence on an outcome variable y. Two predictors will interact in their
#' effect on y, depending on .
#'
#' @param nObs an integer, number of individuals
#' @param nPredictors an integer, number of predictors (i.e., markers/SNPs)
#' @param predType one of c("haploid", "diploid", "continuous", "mixed"),
#' indicating the type of predictor variables to simulate.
#' "haploid": predictors are simulated haploid genotype markers,
#' either 0 or 1, corresponding to the two possible alleles.
#' With LD between neighboring markers according to recombprob.
#' "diploid": predictors are simulated diploid genotype markers. The predictors
#' take values 0, 1, or two, corresponding to homozygotic for the reference,
#' heterozygotic, or homozygotic for the alternate allele. The alternate alleles
#' are sampled according to minMAF. All markers are independent (no LD)
#' "continuous": predictors are normally distributed random variables with
#' means \code{means} and standard deviations \code{SDs}.
#' "mixed": predictors contain both categorical (0/1) and continuous variables,
#' according to \code{propBinary}. Categorical predictors are sampled
#' from 0 or 1 according to \code{minMAF}, and continous are created as
#' described for \code{predType = "continuous"}
#' @param recombprob numeric vector containing values between 0 and 1,
#' representing probablities of recombination; at each successive locus,
#' one of these values is sampled randomly to determine how many of the
#' individuals will "flip" their genotype; the more '0' values in this vector,
#' the more LD there will be between the markers. Higher values will result
#' in more recombination. Defaults to c(0,0,0,0,0.01,0.3). Only used when
#' \code{predType = "haploid"}, ignored otherwise.
#' @param minMAF numeric between 0 and 1, which indicates the minimum allowed
#' minor allele frequency. Defaults to 0.3. Only used for diploid or mixed
#' predictors, ignored otherwise.
#' @param means Vector of desired means for the predictor variables. Only used
#' for continuous or mixed predictors, ignored otherwise.
#' @param SDs  Vector of desired standard deviations for the predictor
#' variables. Only used for continuous or mixed predictors, ignored otherwise.
#' @param propBinary for mixed predictors, a value between 0 and 1, indicating
#' the proportion of predictors that should be binary as opposed to continuous.
#' Defaults to 0.5.
#' @param outcome One of c("binomial", "numeric").
#' "binomial": Outcome variable will be binary 0/1 labels.
#' "numeric": Outcome variable will be continuous values with mean=0 and SD=1
#' @param intType one of c("none", "pure", "modifyer1", "modifyer2","redundant",
#'  "XOR", "synergistic"), indicating type of interaction to be simulated.
#'  Two random predictors x1 and x2 are selected and their interaction effect on
#'  the outcome y modeled according to intType. In case of mixed predictors,
#'  x1 and x2 are selected so that x1 is a binary, and x2 is a continuous
#'  predictor. Effect sizes (ß1, ß2, and ßi)are set with \code{betaMarginals}.
#'  "none": no interaction. The two markers only have marginal effects
#'  (formula: y=ß1\*x1+ß2\*x2)
#'  "pure": interaction only. y=ßi\*x1\*x2
#'  "modifyer1": marginal effect for x1 and interaction (y=ß1\*x1+ßi\*x1\*x2).
#'  "modifyer2": marginal effect for x2 and interaction (y=ß1\*x2+ßi\*x1\*x2).
#'  Essentially the same as modifyer1, except for mixed predictors.
#'  "redundant": the interaction effect negates individual marginal
#'   effects (y=ß1\*x1+ß2\*x2-ßi\*x1\*x2)
#'  "XOR": extremer version of "redundant". Marginal effects are completely
#'  masked by interaction(y=ß1\*x1+ß2\*x2-ßi\*2\*x1\*x2)
#'  "synergistic": both markers have marginal effects as well as interaction
#'  (y=ß1\*x1+ß2\*x2+ßi\*x1\*x2)
#' @param betaInt A numeric vector of length three, indicating marginal and
#' interaction effect sizes of the interacting predictors, see \code{intType}.
#' Defaults to c(1,1,1).
#' @param nMarginals integer, indicating how many predictors should have
#'  marginal effects (in addition to the two interacting predictors).
#' @param betaMarginals numeric vector of length \code{nMarginals},
#'  indicating effect sizes of the marginal effects.
#' @param noiseSD standard deviation of normally distributed noise (mean 0) that
#' is added to the continuous predictor.
#' @return A list with the following components:
#' x a dataframe containing the predictors
#' y the outcome variable
#' intIDS the indexes of the two predictors that interact
#' marginalIDs the indexes of the predictors with additional marginal effects
#' intType same as input intType.
#' predType same as input predType
#' @export
#'
#' @examples
simData = function(nObs,
                   nPredictors,
                   nMarginals,
                   predType = NA,
                   outcome = NA,
                   intType = NA,
                   recombprob=c(0,0,0,0,0.01,0.3), # used for haploid
                   minMAF=0.3, # used for diploid and mixed
                   means=0,
                   SDs = 1, # used for continuous predictors
                   propBinary = 0.5, # used for mixed
                   betaInt = c(1,1,1),
                   betaMarginals = rep(1,nMarginals),
                   noiseSD = 1){
  # add some argument checks here.
  predType = match.arg(predType, c("haploid", "diploid", "continuous", "mixed"))
  outcome = match.arg(outcome, c("binomial", "numeric"))
  intType = match.arg(intType, c("none", "pure", "modifyer1", "modifyer2",
                                 "redundant", "XOR", "synergistic"))
  stopifnot("nObs must be a positive integer" =
              is.numeric(nObs) & nObs>0 & length(nObs) == 1,
            "nPredictors must be a positive integer" =
              is.numeric(nPredictors) & nPredictors>0 & length(nPredictors) == 1,
            "betaInt must be a numeric vector of length 3" =
              length(betaInt) == 3 & is.numeric(betaInt),
            "nMarginals must be a non-negative integer" =
              is.numeric(nMarginals) & nMarginals>=0 & length(nMarginals) == 1,
            "betaMarginals must be a numeric vector of length 1 or nMarginals" =
              is.numeric(betaMarginals) &  length(betaMarginals) == nMarginals,
            "noiseSD must be a numeric" =
              length(noiseSD) == 1 & is.numeric(noiseSD))

  if(predType == "haploid"){
    stopifnot("recombprob must be a numberic vector with values between 0 and 1" =
                is.numeric(recombprob) & all(recombprob <= 1) & all(recombprob >= 0))
  } else if(predType == "diploid"){
    stopifnot("minMAF must be a number between 0 and 1" =
      length(minMAF) == 1 & is.numeric(minMAF) & minMAF <= 1 & minMAF >= 0)

  } else if (predType == "mixed"){
    stopifnot("minMAF must be a number between 0 and 1" =
                length(minMAF) == 1 & is.numeric(minMAF) & minMAF <= 1 & minMAF >= 0,
              "propBinary must be a number between 0 and 1" =
                length(propBinary) == 1 & is.numeric(propBinary) & propBinary <= 1 & propBinary >= 0,
              "means must be a numeric of length 1 or the number of numeric predictors" =
                is.numeric(means) & length(means) %in% c(1,nPredictors*propBinary),
              "SDs must be a numeric of length 1 or the number of numeric predictors" =
                is.numeric(SDs) & length(SDs) %in% c(1,nPredictors*propBinary))
  } else if (predType == "continuous"){
    stopifnot("means must be a numeric of length 1 or the number of numeric predictors" =
      is.numeric(means) & length(means) %in% c(1, nPredictors),
    "SDs must be a numeric of length 1 or the number of numeric predictors" =
      is.numeric(SDs) & length(SDs) %in% c(1, nPredictors))
  }
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
#' @inheritParams simData
#' @return A dataframe with (0/1) genotypes.
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
#' @inheritParams simData
#' @return A dataframe with (0/1/2) genotypes
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
#' @inheritParams simData
#' @return A dataframe with numeric predictors
#'
#' @examples
simContinuous = function(nObs, nPredictors, means, SDs){
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
#' @inheritParams simData
#' @return A dataframe with \code{nBinary} numeric predictors and
#' \code{nNumeric} numeric predictors
#'
#' @examples
simMixed = function(nObs, nBinary, nNumeric,  means, SDs,
                    minMAF){
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


#' @title create numeric Predictor matrix
#'
#' @description
#' A function to simulate an outcome variable with interactions
#' and marginal effects.
#'
#'
#' @inheritParams simData
#' @return a list with the following components:
#' outcome: outcome variable.
#' intIDS: indices of predictors with were randomly selected to interact
#' marginalIDS: indeices of predictors which were randomly selected to
#' have marginal effects.
#'
#' @examples
simulateY = function(X, outcome,
                     predType,
                     intType,
                     betaInt,
                     nMarginals,
                     betaMarginals,
                     noiseSD){
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
