# generate Data
source("R/simulateData.R")
library(randomForest)

# hap Pred, num Outcome
d = simData(nObs = 100, nPredictors = 200, nMarginals = 0,
            predType = "haploid", outcome = "num", intType = "syn")
x = apply(d$x,2,as.numeric)
y = d$y
rf = randomForest::randomForest(x = x, y = y)
save(rf, file = "data/testRF.RData")
tree1 = getTree(rf)

d = simData(nObs = 100, nPredictors = 200, nMarginals = 0,
            predType = "cont", outcome = "bin", intType = "syn")
x = apply(d$x,2,as.numeric)
y = factor(d$y)
rf = randomForest::randomForest(x = x, y = y,)
save(rf, file = "data/testRF2.RData")
tree2 = getTree(rf)
