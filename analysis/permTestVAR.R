# Temporal permutation testing for network-level VAR models.
#__________________________________________________________________________
#
# Author:
#   Tyler Santander (t.santander@psych.ucsb.edu)
#   Institute for Collaborative Biotechnologies
#   Department of Psychological & Brain Sciences
#   University of California, Santa Barbara
#   September 2019
#__________________________________________________________________________

permTestVAR <- function(modelData, modelObject, rngSeed) {
  
  require(ggplot2)
  require(vars)
  
  # Set seed for random number generator.
  
  set.seed(rngSeed)
  
  # Preallocate matrices to hold permuted coefficients.

  randBrainEst     <- matrix(NA, nrow = 5, ncol = 10000)
  randHormoneEst   <- matrix(NA, nrow = 5, ncol = 10000)
  
  randBrainTstat   <- matrix(NA, nrow = 5, ncol = 10000)
  randHormoneTstat <- matrix(NA, nrow = 5, ncol = 10000)
  
  randBrainRsq     <- matrix(NA, nrow = 1, ncol = 10000)
  randHormoneRsq   <- matrix(NA, nrow = 1, ncol = 10000)
  
  randBrainRMSE    <- matrix(NA, nrow = 1, ncol = 10000)
  randHormoneRMSE  <- matrix(NA, nrow = 1, ncol = 10000)
  
  # Loop over permutations, store coefficients.

  for (iPerm in 1:10000) {
  
    randDat           <- modelData
    randDat[,1]       <- sample(modelData[,1])
    randDat[,2]       <- sample(modelData[,2])
    colnames(randDat) <- c('Brain', 'Hormone')
  
    randVAR <- VAR(randDat, p = 2)
    error   <- residuals(randVAR)
    
    randBrainEst[,iPerm]     <- as.numeric(randVAR$varresult[[1]]$coefficients)
    randHormoneEst[,iPerm]   <- as.numeric(randVAR$varresult[[2]]$coefficients)
    
    randBrainTstat[,iPerm]   <- as.numeric(summary(randVAR)$varresult[[1]]$coefficients[, 't value'])
    randHormoneTstat[,iPerm] <- as.numeric(summary(randVAR)$varresult[[2]]$coefficients[, 't value'])
    
    randBrainRsq[iPerm]      <- summary(randVAR)$varresult[[1]]$r.squared
    randHormoneRsq[iPerm]    <- summary(randVAR)$varresult[[2]]$r.squared
    
    randBrainRMSE[iPerm]     <- sqrt((t(error[,1]) %*% error[,1]) / length(error[,1]))
    randHormoneRMSE[iPerm]   <- sqrt((t(error[,2]) %*% error[,2]) / length(error[,2]))
  
  }
  
  # Get two-tailed empirical p-values.

  trueBrainEst      <- as.numeric(modelObject$varresult[[1]]$coefficients)
  trueHormoneEst    <- as.numeric(modelObject$varresult[[2]]$coefficients)
  
  trueBrainTstat    <- as.numeric(summary(modelObject)$varresult[[1]]$coefficients[, 't value'])
  trueHormoneTstat  <- as.numeric(summary(modelObject)$varresult[[2]]$coefficients[, 't value'])
  
  trueBrainRsq      <- summary(modelObject)$varresult[[1]]$r.squared
  trueHormoneRsq    <- summary(modelObject)$varresult[[2]]$r.squared
  
  error             <- residuals(modelObject)
  trueBrainRMSE     <- sqrt((t(error[,1]) %*% error[,1]) / length(error[,1]))
  trueHormoneRMSE   <- sqrt((t(error[,2]) %*% error[,2]) / length(error[,2]))
  
  pBrainEst         <- matrix(NA, nrow = 5, ncol = 1)
  pHormoneEst       <- matrix(NA, nrow = 5, ncol = 1)
  
  pBrainTstat       <- matrix(NA, nrow = 5, ncol = 1)
  pHormoneTstat     <- matrix(NA, nrow = 5, ncol = 1)
  
  for (iCoef in 1:5) {
    
    pBrainEst[iCoef]     <- sum(abs(randBrainEst[iCoef,]) >= abs(trueBrainEst[iCoef])) / 10000
    pHormoneEst[iCoef]   <- sum(abs(randHormoneEst[iCoef,]) >= abs(trueHormoneEst[iCoef])) / 10000
    
    pBrainTstat[iCoef]   <- sum(abs(randBrainTstat[iCoef,]) >= abs(trueBrainTstat[iCoef])) / 10000
    pHormoneTstat[iCoef] <- sum(abs(randHormoneTstat[iCoef,]) >= abs(trueHormoneTstat[iCoef])) / 10000
    
  }
  
  pBrainRsq    <- sum(randBrainRsq >= trueBrainRsq) / 10000
  pHormoneRsq  <- sum(randHormoneRsq >= trueHormoneRsq) / 10000
  
  pBrainRMSE   <- sum(randBrainRMSE <= trueBrainRMSE[1]) / 10000
  pHormoneRMSE <- sum(randHormoneRMSE <= trueHormoneRMSE[1]) / 10000
  
  # Plot null distributions for parameters on brain outcome.
  
  randBrainEst           <- data.frame(t(randBrainEst))
  colnames(randBrainEst) <- names(modelObject$varresult[[1]]$coefficients)
  
  for (iCoef in 1:5) {
    
    print(ggplot(randBrainEst, aes(x=randBrainEst[,iCoef])) + 
            geom_histogram(aes(y=..density..),
                           binwidth=.05,
                           colour="black", fill="white") +
            geom_density(alpha=.2, fill="skyblue") +
            geom_vline(aes(xintercept=trueBrainEst[iCoef]),
                       color="red", linetype="dashed", size=1) +
            xlab(paste(names(randVAR$varresult[1]), ': ', colnames(randBrainEst)[iCoef], ' (p = ', pBrainEst[iCoef], ')', sep = '')) + 
            ylab('Density') +
            theme_Publication())
    
  }
  
  # Plot null distributions for parameters on hormone outcome.
  
  randHormoneEst           <- data.frame(t(randHormoneEst))
  colnames(randHormoneEst) <- names(modelObject$varresult[[2]]$coefficients)
  
  for (iCoef in 1:5) {
    
    print(ggplot(randHormoneEst, aes(x=randHormoneEst[,iCoef])) + 
            geom_histogram(aes(y=..density..),
                           binwidth=.05,
                           colour="black", fill="white") +
            geom_density(alpha=.2, fill="skyblue") +
            geom_vline(aes(xintercept=trueHormoneEst[iCoef]),
                       color="red", linetype="dashed", size=1) +
            xlab(paste(names(randVAR$varresult[2]), ': ', colnames(randHormoneEst)[iCoef], ' (p = ', pHormoneEst[iCoef], ')', sep = '')) + 
            ylab('Density') +
            theme_Publication())
    
  }
  
  # Plot null distributions for test statistics on brain outcome.
  
  randBrainTstat           <- data.frame(t(randBrainTstat))
  colnames(randBrainTstat) <- names(modelObject$varresult[[1]]$coefficients)
  
  for (iCoef in 1:5) {
    
    print(ggplot(randBrainTstat, aes(x=randBrainTstat[,iCoef])) + 
            geom_histogram(aes(y=..density..),
                           binwidth=.05,
                           colour="black", fill="white") +
            geom_density(alpha=.2, fill="skyblue") +
            geom_vline(aes(xintercept=trueBrainTstat[iCoef]),
                       color="red", linetype="dashed", size=1) +
            xlab(paste(names(randVAR$varresult[1]), ': ', colnames(randBrainTstat)[iCoef], ' (p = ', pBrainTstat[iCoef], ')', sep = '')) + 
            ylab('Density') +
            theme_Publication())
    
  }
  
  # Plot null distributions for test statistics on hormone outcome.
  
  randHormoneTstat           <- data.frame(t(randHormoneTstat))
  colnames(randHormoneTstat) <- names(modelObject$varresult[[2]]$coefficients)
  
  for (iCoef in 1:5) {
    
    print(ggplot(randHormoneTstat, aes(x=randHormoneTstat[,iCoef])) + 
            geom_histogram(aes(y=..density..),
                           binwidth=.05,
                           colour="black", fill="white") +
            geom_density(alpha=.2, fill="skyblue") +
            geom_vline(aes(xintercept=trueHormoneTstat[iCoef]),
                       color="red", linetype="dashed", size=1) +
            xlab(paste(names(randVAR$varresult[2]), ': ', colnames(randHormoneTstat)[iCoef], ' (p = ', pHormoneTstat[iCoef], ')', sep = '')) + 
            ylab('Density') +
            theme_Publication())
    
  }
  
  # Plot null distributions for Rsquared.
  
  randBrainRsq <- data.frame(t(randBrainRsq))
  
  print(ggplot(randBrainRsq, aes(x=randBrainRsq[,1])) + 
          geom_histogram(aes(y=..density..),
                         binwidth=.05,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="skyblue") +
          geom_vline(aes(xintercept=trueBrainRsq),
                     color="red", linetype="dashed", size=1) +
          xlab(paste(names(randVAR$varresult[1]), ': Rsquared (p = ', pBrainRsq, ')', sep = '')) + 
          ylab('Density') +
          theme_Publication())
  
  randHormoneRsq <- data.frame(t(randHormoneRsq))
  
  print(ggplot(randHormoneRsq, aes(x=randHormoneRsq[,1])) + 
          geom_histogram(aes(y=..density..),
                         binwidth=.05,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="skyblue") +
          geom_vline(aes(xintercept=trueHormoneRsq),
                     color="red", linetype="dashed", size=1) +
          xlab(paste(names(randVAR$varresult[2]), ': Rsquared (p = ', pHormoneRsq, ')', sep = '')) + 
          ylab('Density') +
          theme_Publication())
  
  # Plot null distributions for RMSE.
  
  randBrainRMSE <- data.frame(t(randBrainRMSE))
  
  print(ggplot(randBrainRMSE, aes(x=randBrainRMSE[,1])) + 
          geom_histogram(aes(y=..density..),
                         binwidth=.05,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="skyblue") +
          geom_vline(aes(xintercept=trueBrainRMSE[1]),
                     color="red", linetype="dashed", size=1) +
          xlab(paste(names(randVAR$varresult[1]), ': RMSE (p = ', pBrainRMSE, ')', sep = '')) + 
          ylab('Density') +
          theme_Publication())
  
  randHormoneRMSE <- data.frame(t(randHormoneRMSE))
  
  print(ggplot(randHormoneRMSE, aes(x=randHormoneRMSE[,1])) + 
          geom_histogram(aes(y=..density..),
                         binwidth=.05,
                         colour="black", fill="white") +
          geom_density(alpha=.2, fill="skyblue") +
          geom_vline(aes(xintercept=trueHormoneRMSE[1]),
                     color="red", linetype="dashed", size=1) +
          xlab(paste(names(randVAR$varresult[2]), ': RMSE (p = ', pHormoneRMSE, ')', sep = '')) + 
          ylab('Density') +
          theme_Publication())
  
  # Collect results into a data frame for storage or printing to the command line.
  
  Model     <- c(rep(names(randVAR$varresult[1]),5), rep(names(randVAR$varresult[2]), 5))
  Term      <- c(names(modelObject$varresult[[1]]$coefficients), names(modelObject$varresult[[2]]$coefficients))
  Estimate  <- c(round(trueBrainEst,3), round(trueHormoneEst,3))
  pValue    <- c(pBrainEst, pHormoneEst)
  permEst   <- data.frame(Model, Term, Estimate, pValue)
  
  Tstat     <- c(round(trueBrainTstat,3), round(trueHormoneTstat,3))
  pValue    <- c(pBrainTstat, pHormoneTstat)
  permTstat <- data.frame(Model, Term, Tstat, pValue)
  
  Model     <- c(names(randVAR$varresult[1]), names(randVAR$varresult[2]))
  Rsq       <- c(round(trueBrainRsq,3), round(trueHormoneRsq,3))
  pValue    <- c(pBrainRsq, pHormoneRsq)
  permRsq   <- data.frame(Model, Rsq, pValue)
  
  RMSE      <- c(round(trueBrainRMSE,3), round(trueHormoneRMSE,3))
  pValue    <- c(pBrainRMSE, pHormoneRMSE)
  permRMSE  <- data.frame(Model, RMSE, pValue)
  
  permResult <- list('Coefficients' = permEst, 'Tstat' = permTstat, 
                     'Rsquared' = permRsq, 'RMSE' = permRMSE)
  
  return(permResult)
  
}

theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.border = element_blank(),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.text.x = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.margin = margin(0.01, 0.01, 0.01, 0.01, unit= "cm"),
            legend.title = element_text(face="italic"),
            plot.margin = unit(c(10,5,5,5),"mm"),
            plot.background = element_blank(),
            strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}