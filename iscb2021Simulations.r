#R code for simulations reported in ISCB 2021 presentation
#Reference based multiple imputation for trials - what's the right variance and how to estimate it?
#Jonathan Bartlett

library(dejaVu)
library(bootImpute)

#perform MI under the copy reference assumption
impM <- function(obsData,M,proper=TRUE) {
  library(dejaVu)
  #first import obsData. Make ids all unique
  obsData$Id <- 1:dim(obsData)[1]
  obsData <- ImportSim(MakeDejaData(subset(obsData,select=c("Id", "arm")), arm="arm", Id="Id"),
                       event.times=expandEventCount(count=obsData$observed.events, time=obsData$censored.time),
                       status="dropout",
                       study.time=365, censored.time=obsData$censored.time, allow.beyond.study=FALSE)
  #fit model to observed data
  fit <- Simfit(obsData)
  #impute using copy reference
  imps <- Impute(fit, copy_reference(proper=proper), M)
  imps_list <- vector(mode = "list", length = M)
  for (i in 1:M) {
    imps_list[[i]] <- GetImputedDataSet(imps,index=i)
  }
  imps_list
}

analyseImp <- function(singleImpData) {
  library(dejaVu)
  fitted <- Simfit(singleImpData,family="negbin")
  crMIEst <- log(summary(fitted)$treatment.effect)
  crMIEst
}
  
#define function to run simulations
runSim <- function(nSim=100,nBoot=500,dropoutRate=0.0025) {

  rubin <- array(0, dim=c(nSim,4))
  vonHippel <- array(0, dim=c(nSim,4))
  
  for (i in 1:nSim) {
    print(i)
    
    complete <- SimulateComplete(study.time=365,
                                 number.subjects=500,
                                 event.rates=c(0.01,0.005),
                                 dispersions=0.25)
    #impose some dropout
    with.MCAR.dropout <- SimulateDropout(complete,
                                         drop.mechanism=ConstantRateDrop(rate=dropoutRate,
                                                                         var=0))
    
    fit <- Simfit(with.MCAR.dropout)
    imps <- Impute(fit, copy_reference(), 10)
    impfitted <- Simfit(imps,
                     family="negbin")
    result <- summary(impfitted)
    rubin[i,] <- c(log(result$treatment.effect),
                   result$se,
                   log(result$treatment.effect)-qt(0.975,df=result$adjusted.df)*result$se,
                   log(result$treatment.effect)+qt(0.975,df=result$adjusted.df)*result$se)
  
    #von Hippel and Bartlett
    
    #save seed because with multiple cores bootImpute requires us to set the seed
    oldseed <- .Random.seed
    runBootImpute <- bootImpute(with.MCAR.dropout$data, impM, nBoot=nBoot, nImp=2, M=2, proper=FALSE,
                                nCores=10, seed=4126)
    #restore seed to where it was
    .Random.seed <- oldseed
  
    runBootImputeAnalyse <- bootImputeAnalyse(runBootImpute, analyseImp,nCores=1, quiet=TRUE)
    vonHippel[i,] <- c(runBootImputeAnalyse$ests,
                       runBootImputeAnalyse$var^0.5,
                       runBootImputeAnalyse$ci)
    
  }
  
  list(rubin=rubin, vonHippel=vonHippel)
}

sumSim <- function(simSet) {
  print("Rubin")
  print(c(mean(simSet$rubin[,1]),
    sd(simSet$rubin[,1]),
    mean(simSet$rubin[,2])))
  print("von Hippel")
  print(c(mean(simSet$vonHippel[,1]),
          sd(simSet$vonHippel[,1]),
          mean(simSet$vonHippel[,2])))
}

#probability of not dropping out before end of fup is
exp(-365*0.00025)
exp(-365*0.0025)

set.seed(72347218)
nSim <- 1000
run1 <- runSim(nSim=nSim, nBoot=200, dropoutRate=0.00025)
run2 <- runSim(nSim=nSim, nBoot=200, dropoutRate=0.0025)

sumSim(run1)
sumSim(run2)

#construct results table for Beamer presentation
resTable <- array(0, dim=c(4,4))
row.names(resTable) <- rep(c("Rubin", "Bootstrap"),2)
colnames(resTable) <- c("Est. log", "Emp.  SE", "Est. SE", "SE ratio")
resTable[1,1:3] <- c(mean(run1$rubin[,1]),
                  sd(run1$rubin[,1]),
                  mean(run1$rubin[,2]))
resTable[2,1:3] <- c(mean(run1$vonHippel[,1]),
                  sd(run1$vonHippel[,1]),
                  mean(run1$vonHippel[,2]))
resTable[3,1:3] <- c(mean(run2$rubin[,1]),
                  sd(run2$rubin[,1]),
                  mean(run2$rubin[,2]))
resTable[4,1:3] <- c(mean(run2$vonHippel[,1]),
                  sd(run2$vonHippel[,1]),
                  mean(run2$vonHippel[,2]))
resTable[,4] <- resTable[,3]/resTable[,2]
resTable

library(xtable)
xtable(resTable, digits=c(2,2,3,3,2))