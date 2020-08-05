# Temporal permutation testing for edgewise VAR models.
#__________________________________________________________________________
#
# Author:
#   Tyler Santander (t.santander@psych.ucsb.edu)
#   Institute for Collaborative Biotechnologies
#   Department of Psychological & Brain Sciences
#   University of California, Santa Barbara
#   September 2019
#__________________________________________________________________________

permTestEdgeVAR <- function(modelData, modelObject, rngSeed) {
  
  require(vars)
  
  # Set seed for random number generator.
  
  set.seed(rngSeed)
  
  # Preallocate matrices to hold permuted coefficients.
  
  randBrainEst    <- matrix(NA, nrow = 5, ncol = 10000)
  randHormoneEst  <- matrix(NA, nrow = 5, ncol = 10000)
  
  randBrainRsq    <- matrix(NA, nrow = 1, ncol = 10000)
  randHormoneRsq  <- matrix(NA, nrow = 1, ncol = 10000)
  
  randBrainRMSE   <- matrix(NA, nrow = 1, ncol = 10000)
  randHormoneRMSE <- matrix(NA, nrow = 1, ncol = 10000)
  
  # Loop over permutations, store coefficients.
  
  for (iPerm in 1:10000) {
    
    randDat     <- modelData
    randDat[,1] <- sample(modelData[,1])
    randDat[,2] <- sample(modelData[,2])
    
    randVAR <- VAR(randDat, p = 2)
    error   <- residuals(randVAR)
    
    randBrainEst[,iPerm]   <- as.numeric(summary(randVAR)$varresult[[1]]$coefficients[, 't value'])
    randHormoneEst[,iPerm] <- as.numeric(summary(randVAR)$varresult[[2]]$coefficients[, 't value'])
    
    randBrainRsq[iPerm]    <- summary(randVAR)$varresult[[1]]$r.squared
    randHormoneRsq[iPerm]  <- summary(randVAR)$varresult[[2]]$r.squared
    
    randBrainRMSE[iPerm]   <- sqrt((t(error[,1]) %*% error[,1]) / length(error[,1]))
    randHormoneRMSE[iPerm] <- sqrt((t(error[,2]) %*% error[,2]) / length(error[,2]))
    
  }
  
  # Get two-tailed empirical p-values.
  
  trueBrainEst    <- as.numeric(summary(modelObject)$varresult[[1]]$coefficients[, 't value'])
  trueHormoneEst  <- as.numeric(summary(modelObject)$varresult[[2]]$coefficients[, 't value'])
  
  trueBrainRsq    <- summary(modelObject)$varresult[[1]]$r.squared
  trueHormoneRsq  <- summary(modelObject)$varresult[[2]]$r.squared
  
  error           <- residuals(modelObject)
  trueBrainRMSE   <- sqrt((t(error[,1]) %*% error[,1]) / length(error[,1]))
  trueHormoneRMSE <- sqrt((t(error[,2]) %*% error[,2]) / length(error[,2]))
  
  pBrainEst   <- matrix(NA, nrow = 5, ncol = 1)
  pHormoneEst <- matrix(NA, nrow = 5, ncol = 1)
  
  for (iCoef in 1:5) {
    
    pBrainEst[iCoef]   <- sum(abs(randBrainEst[iCoef,]) >= abs(trueBrainEst[iCoef])) / 10000
    pHormoneEst[iCoef] <- sum(abs(randHormoneEst[iCoef,]) >= abs(trueHormoneEst[iCoef])) / 10000
    
  }
  
  pBrainRsq    <- sum(randBrainRsq >= trueBrainRsq) / 10000
  pHormoneRsq  <- sum(randHormoneRsq >= trueHormoneRsq) / 10000
  
  pBrainRMSE   <- sum(randBrainRMSE <= trueBrainRMSE[1]) / 10000
  pHormoneRMSE <- sum(randHormoneRMSE <= trueHormoneRMSE[1]) / 10000
  
  # Collect results into a data frame for storage or printing to the command line.
  
  Model    <- c(rep(names(randVAR$varresult[1]),5), rep(names(randVAR$varresult[2]), 5))
  Term     <- c(names(modelObject$varresult[[1]]$coefficients), names(modelObject$varresult[[2]]$coefficients))
  Estimate <- c(round(trueBrainEst,3), round(trueHormoneEst,3))
  pValue   <- c(pBrainEst, pHormoneEst)
  permEst  <- data.frame(Model, Term, Estimate, pValue)
  
  Model    <- c(names(randVAR$varresult[1]), names(randVAR$varresult[2]))
  Rsq      <- c(round(trueBrainRsq,3), round(trueHormoneRsq,3))
  pValue   <- c(pBrainRsq, pHormoneRsq)
  permRsq  <- data.frame(Model, Rsq, pValue)
  
  RMSE     <- c(round(trueBrainRMSE[1],3), round(trueHormoneRMSE[1],3))
  pValue   <- c(pBrainRMSE, pHormoneRMSE)
  permRMSE <- data.frame(Model, RMSE, pValue)
  
  permResult <- list('Coefficients' = permEst, 'Rsquared' = permRsq, 'RMSE' = permRMSE)
  
  return(permResult)
  
}