source("R/simulateData.R")
library(randomForest)

# hap Pred, num Outcome
d = simData(nObs = 100, nPredictors = 200, nMarginals = 0,
            predType = "haploid", outcome = "num", intType = "syn")
x = apply(d$x,2,as.numeric)
y = d$y
rf1 = randomForest::randomForest(x = x, y = y, keep.inbag = T)
rf2 = RandomForestExtended::randomForest(x = x, y = y, keep.inbag = T)

compareRootPreds = function(rfObj){
  t(sapply(1:rfObj$ntree, function(i){
    pred = rfObj$forest$nodepred[1,i]
    oobPred = rfObj$forest$oobpred[1,i]
    ib = rfObj$inbag[,i]
    ibData = rep(1:length(y),ib)
    ibY = mean(y[ibData])
    obY = mean(y[ib==0])
    return(c(pred = pred,oobPred = oobPred, ibY= ibY, obY=obY))
  }))
}
compareNode2Preds = function(rfObj){
  t(sapply(1:rfObj$ntree, function(i){
    leftN = rfObj$forest$leftDaughter[1,i]
    pred = rfObj$forest$nodepred[leftN,i]
    oobPred = rfObj$forest$oobpred[leftN,i]
    splitV = rfObj$forest$bestvar[1,i]
    splitP = rfObj$forest$xbestsplit[1,i]

    ib = rfObj$inbag[,i]
    ibData = rep(1:length(y),ib)
    ibY = mean(y[ibData][x[ibData,splitV] < splitP])
    obY = mean(y[ib==0][x[ib==0,splitV] < splitP])
    return(c(pred = pred,oobPred = oobPred, ibY= ibY, obY=obY))
  }))
}
res1 = compareRootPreds(rfObj = rf1)
res2 = compareRootPreds(rfObj = rf2)
plot(as.data.frame(res1), main = "randomForest")
plot(as.data.frame(res2), main = "RandomForestExtended")


res3 = compareNode2Preds(rfObj = rf1)
res4 = compareNode2Preds(rfObj = rf2)
plot(as.data.frame(res3), main = "randomForest")
plot(as.data.frame(res4), main = "RandomForestExtended")

