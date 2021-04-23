library(ggplot2)
library(gridExtra)
source("R/simulateData.R")
source("R/helperFunctions.R")

nObs =1000
nPredictors =1000
# betaMarginals = rep(0.5,nMarginals)
intTypes = setNames(c("none (y=x1+x2)", "pure (y=x1*x2)", "modifyer1 (y=x1+x1*x2)",
                      "modifyer2 (y=x2+x1*x2)",
                      "redundant (y=x1+x2-x1*x2)", "XOR (y=x1+x2-2*x1*x2)",
                      "synergistic (y=x1+x2+x1*x2)"),
                    c("none", "pure", "modifyer1", "modifyer2",
                      "redundant", "XOR", "synergistic"))
outTypes = c("numeric", "binomial")
predTypes =  c("haploid", "diploid", "continuous", "mixed")

# haploid predictors, numeric outcome #####
pdfandpng(path = "fig/simData_haplPred_numOutcome",
          height = 9, width = 6,
          expr = expr({
            par(mfrow = c(4,2), mar = c(3,3,2,0.5), mgp = c(2.1,1,0))
            for(i in names(intTypes)){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "haploid",
                          outcome = "numeric",
                          intType = i,
                          nMarginals = 0, noiseSD = 0.2)
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])
              boxplot(y ~ x1:x2, data = d,  las = 1,
                      main = intTypes[i], ylab = "y")
              mod = lm(y ~ x1*x2, data = d)
              p = summary(mod)$coefficients [2:4,4]
              legend("topleft", bty = "n", cex = 0.8,
                     legend = paste0(c("p(x1) = ", "p(x2) = ", "p(x1:x2) = "),
                                     format(p, digits = 2)))
            }
          }))
#####

# diploid predictors, numeric outcome #####
pdfandpng(path = "fig/simData_diplPred_numOutcome",
          height = 9, width = 6,
          expr = expr({
            par(mfrow = c(4,2), mar = c(3,3,1,0.5))
            for(i in names(intTypes)){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "diploid",
                          outcome = "numeric",
                          intType = i,
                          nMarginals = 0, noiseSD = 0.5)
              d = data.frame(y = D$y,
                             x1 = D$x[,D$intIDs[1]],
                             x2 = D$x[,D$intIDs[2]])
              boxplot(y ~ x1:x2, data = d,  las = 1, mgp = c(2,0.7,0),
                      main = intTypes[i], names = NA,
                      col = c("darkgoldenrod2", "chocolate3", "darkred"))
              abline(v = c(3.5, 6.5))
              mtext(text = c("bb", "bB", "BB"), side = 1, at = c(2,5,8),
                    line= 2, col = c("deepskyblue", "royalblue3", "blue"))
              arrows(x0 = c(0.55, 3.55, 6.55), x1 =  c(3.45, 6.45, 9.45),
                     y0 = (min(d$y)-(max(d$y)-min(d$y))*0.2),length = 0.05,
                     col = c("deepskyblue", "royalblue3", "blue"),
                     angle = 90,  code = 3, xpd = NA)
              mtext(text = c("aa", "aA", "AA"), side = 1,
                    at = 1:9, line= 0.8,
                    col = c("darkgoldenrod2", "chocolate3", "darkred"))
              mod = lm(y ~ x1*x2, data = d)
              p = summary(mod)$coefficients [2:4,4]
              legend("topleft", bty = "n",
                     legend = paste0(c("p(x1) = ", "p(x2) = ", "p(x1:x2) = "),
                                     format(p, digits = 2)))
            }
            }))
#####

# numeric predictors, numeric outcome #####
pdfandpng(path = "fig/simData_numPred_numOutcome",
          height = 9, width = 6,
          expr = expr({
            plotList = lapply(names(intTypes), function(i){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "continuous",
                          outcome = "numeric", nMarginals = 0,
                          intType = i, betaInt = c(10,10,50))
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])
              d = d[abs(d$x1) < 1.5 & abs(d$x2) < 1.5,]
              g = ggplot(d) +
                geom_point(mapping = aes(x = x1, y = x2, col = y), size = 2)+
                ggtitle(intTypes[i]) +
                labs(col = "y")
              return(g)
            })
            grid.arrange(grobs = plotList, ncol = 2)
          }))

#####

# mixed predictors, numeric outcome #####
pdfandpng(path = "fig/simData_mixPred_numOutcome",
          height = 9, width = 6,
          expr = expr({
            par(mfrow = c(4,2), mar = c(3,3,1.5,0.5))
            for(i in names(intTypes)){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "mixed",
                          outcome = "numeric",
                          intType = i,
                          nMarginals = 0, noiseSD = 0.2)
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])
              plot(y~x2,data = d, col = d$x1+1, main = intTypes[i],
                   xlab = "x2", ylab = "y", las = 1, mgp = c(2,0.7,0))
              mod = lm(y~x1*x2, data = d)
              p = summary(mod)$coefficients [2:4,4]
              legend("topleft", bty = "n",
                     legend = paste0(c("p(x1) = ", "p(x2) = ", "p(x1:x2) = "),
                                     format(p, digits = 2)))
              legend("right", col = 1:2, legend = 0:1, pch = 1, title = "x1")
            }
          }))

#####

# haploid predictors, categorial outcome #####
pdfandpng(path = "fig/simData_haplPred_catOutcome",
          height = 9, width = 6,
          expr = expr({
            plotList = lapply(names(intTypes), function(i){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "haploid",
                          outcome = "binomial",
                          intType = i,
                          nMarginals = 0, noiseSD = 0)
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])
              temp = split(d$y, interaction(d$x1, d$x2))
              pD = data.frame(do.call(rbind,strsplit(names(temp), split = "[.]")),
                              prop = sapply(temp, mean))
              g = ggplot(data = pD, mapping = aes(x = X1, y = X2))+
                geom_raster(aes(fill = prop)) +
                ggtitle(intTypes[i])
              return(g)
            })
            grid.arrange(grobs = plotList, ncol = 2)
          }))

######

# diploid predictors, categorial outcome #####
pdfandpng(path = "fig/simData_diplPred_catOutcome",
          height = 9, width = 6,
          expr = expr({
            plotList = lapply(names(intTypes), function(i){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "diploid",
                          outcome = "binomial",
                          intType = i,
                          nMarginals = 0, noiseSD = 0)
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])
              d$x1 = c("aa", "aA", "AA")[d$x1+1]
              d$x2 = c("bb", "bB", "BB")[d$x2+1]
              temp = split(d$y, interaction(d$x1, d$x2))
              pD = data.frame(do.call(rbind,strsplit(names(temp), split = "[.]")),
                              prop = sapply(temp, mean))
              g = ggplot(data = pD, mapping = aes(x = X1, y = X2))+
                geom_raster(aes(fill = prop)) +
                ggtitle(intTypes[i])
              return(g)
            })
            grid.arrange(grobs = plotList, ncol = 2)
          }))

######

# num predictors, categorical outcome ######
pdfandpng(path = "fig/simData_numPred_catOutcome",
          height = 9, width = 6,
          expr = expr({
            colFunc = rgb(red = c(0,1), green = c(0,0),blue =
                            c(0,0), c(0.3,0.3), maxColorValue = 1)
            par(mfrow = c(4,2), mar = c(3,3,2,0.5))
            for(i in names(intTypes)){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "continuous",
                          outcome = "binomial",
                          intType = i,betaInt = c(10,10,10),
                          nMarginals = 0, noiseSD = 0)
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])

              plot(d$x1, d$x2, col = colFunc[d$y+1], pch = 19, cex = 1,
                   las = 1, mgp = c(2,0.7,0), main = intTypes[i])
              mod = glm(y ~ x1*x2, d, family=binomial(link='logit'))
              p = summary(mod)$coefficients [2:4,4]
              legend("topleft",
                     legend = paste0(c("p(x1) = ", "p(x2) = ", "p(x1:x2) = "),
                                     format(p, digits = 2)))
            }
          }))

######

# mixed predictors, categorial outcome #####
pdfandpng(path = "fig/simData_mixPred_catOutcome",
          height = 9, width = 6,
          expr = expr({
            par(mfrow = c(4,2), mar = c(3,3,1,0.5))
            for(i in names(intTypes)){
              D = simData(nObs = nObs, nPredictors = nPredictors,
                          predType = "mixed",
                          outcome = "binomial",
                          intType = i,
                          nMarginals = 0, noiseSD = 0)
              d = data.frame(y = D$y, x1 = D$x[,D$intIDs[1]], x2 = D$x[,D$intIDs[2]])
              plot(d$y~d$x2, col = d$x1+1, main = intTypes[i],
                   xlab = "x2", ylab = "y", las = 1, mgp = c(2,0.7,0))
              subD0 = d[d$x1 == 0,]
              subD0 = subD0[order(subD0$x2),]
              mod0 = glm(y~x2, subD0, family = "binomial")
              lines(subD0$x2,mod0$fitted.values, col ="black")
              subD1 = d[d$x1 == 1,]
              subD1 = subD1[order(subD1$x2),]
              mod1 = glm(y~x2,subD1, family = "binomial")
              lines(subD1$x2,mod1$fitted.values, col ="red")
              legend("left", col = c("black", "red"), legend = 0:1,
                     lty = 1, pch = 1,title = "x1")
            }
          }))

#####

