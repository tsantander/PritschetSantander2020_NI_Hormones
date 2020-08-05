# Estimate second-order vector autoregressive models of edgwise coherence 
# vs. estradiol using parallel processing.
#__________________________________________________________________________
#
# Author:
#   Tyler Santander (t.santander@psych.ucsb.edu)
#   Institute for Collaborative Biotechnologies
#   Department of Psychological & Brain Sciences
#   University of California, Santa Barbara
#   September 2019
#__________________________________________________________________________

# Load in required packages.
#--------------------------------------------------------------------------

require(vars)
require(doMC)

# Specify number of cores for parallel estimation.
#--------------------------------------------------------------------------

registerDoMC(cores=24)

# Navigate to where the data live, load things into memory.
#--------------------------------------------------------------------------

setwd('/home/tyler/28andMe/results.network.rest') 
source('./permTestEdgeVAR.R')

zCoh <- read.csv('./edgewiseCoherence.csv', header = FALSE)
dat  <- read.csv('./data.efficiency.csv', header = TRUE)

zCoh           <- data.frame(scale(zCoh))
zDat           <- data.frame(scale(dat[,2:11]))
colnames(zDat) <- colnames(dat)[2:11]

nEdge <- dim(zCoh)[2]

# Loop over edges in parallel, fit models.
#--------------------------------------------------------------------------

nptOut <- foreach(iEdge=1:nEdge, .combine=rbind, .packages='vars') %dopar% {
  
  modelDat           <- data.frame(zCoh[,iEdge], zDat$Estro)
  colnames(modelDat) <- c('Edge', 'Estro')
  
  edgeVAR <- VAR(modelDat, p = 2)
  
  permResult <- permTestEdgeVAR(modelDat, edgeVAR, sample(1:5e5, 1))
  paramOut   <- permResult$Coefficients
  
  edgeCoef  <- paramOut$Estimate[paramOut$Model == 'Edge']
  estroCoef <- paramOut$Estimate[paramOut$Model == 'Estro']
  edgeNPT   <- paramOut$pValue[paramOut$Model == 'Edge']
  estroNPT  <- paramOut$pValue[paramOut$Model == 'Estro']
  
  npt <- c(edgeCoef, edgeNPT, estroCoef, estroNPT)
  
}

# Save output.
#--------------------------------------------------------------------------

save(nptOut, file = 'estroEdgeVAR.rda', compress = 'xz')
