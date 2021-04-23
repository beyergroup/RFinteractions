source("R/simulateData.R")


intTypes = c("none", "pure", "modifyer1", "modifyer2",
             "redundant", "XOR", "synergistic")
outTypes = c("numeric", "binomial")
predTypes =  c("haploid", "diploid", "continuous", "mixed")
combinations = expand.grid(predTypes, outTypes, intTypes, stringsAsFactors = F)
colnames(combinations) = c("predType", "outType", "intType")
nObs =200
nPredictors =1000
nSimulations = 10
nMarginals = 5
dir.create("data/simulations")
set.seed(34512)
dumpVar = apply(combinations,1, function(x){
  dir.create(paste0("data/simulations/",
                    "pred_",x["predType"],
                    "_out_", x["outType"],
                    "_int_", x["intType"]))
  temp = sapply(1:nSimulations, function(i){
    d = simData(nObs = nObs,
                nPredictors = nPredictors,
                propBinary = 0.5, # proportion of binary for mixed predictors
                means=0, SDs = 1, # mean and sd for continuous predictors
                recombprob=c(0,0,0,0,0.01,0.3), # used for haploid
                minMAF=0.3, # used for diploid and mixed
                predType = x["predType"],
                outcome = x["outType"],
                intType = x["intType"], betaInt = c(1,1,1),
                nMarginals = nMarginals, betaMarginals = rep(1,nMarginals),
                noiseSD = 0.5)
    save(d, file = paste0("data/simulations/",
                          "pred_",x["predType"],
                          "_out_", x["outType"],
                          "_int_", x["intType"],
                          "/sim",i,".RData"))
  })
})
