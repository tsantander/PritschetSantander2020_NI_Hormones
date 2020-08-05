# Estimate second-order vector autoregressive models for within-network 
# efficiency and between-network participation vs. estradiol. 
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

library(tidyverse)
library(vars)

# Navigate to where the data live, load some additional functions.
#--------------------------------------------------------------------------

setwd('~/Documents/28andMe/results.network.rest')
source('./permTestVAR.R')

# Fit models using within-network efficiency.
#--------------------------------------------------------------------------

# Load in the data.

dat <- read.csv('./data.efficiency.csv')

# Z-transform.

zDat           <- data.frame(scale(dat[,2:27]))
colnames(zDat) <- colnames(dat)[2:27]

# Run VAR for FCN:

fcnDat           <- data.frame(zDat$Control, zDat$Estro)
colnames(fcnDat) <- c('FCN', 'Estradiol')

fcnVAR <- VAR(fcnDat, p = 2)
summary(fcnVAR)

fcnPerm <- permTestVAR(fcnDat, fcnVAR, 932)

# Run VAR for DMN:

dmnDat           <- data.frame(zDat$DMN, zDat$Estro)
colnames(dmnDat) <- c('DMN', 'Estradiol')

dmnVAR <- VAR(dmnDat, p = 2)
summary(dmnVAR)

dmnPerm <- permTestVAR(dmnDat, dmnVAR, 235)

# Run VAR for DAN:

danDat           <- data.frame(zDat$DorsAttn, zDat$Estro)
colnames(danDat) <- c('DorsAttn', 'Estradiol')

danVAR <- VAR(danDat, p = 2)
summary(danVAR)

danPerm <- permTestVAR(danDat, danVAR, 263)

# Run VAR for Limbic:

limbicDat           <- data.frame(zDat$Limbic, zDat$Estro)
colnames(limbicDat) <- c('Limbic', 'Estradiol')

limbicVAR <- VAR(limbicDat, p = 2)
summary(limbicVAR)

limbicPerm <- permTestVAR(limbicDat, limbicVAR, 134)

# Run VAR for VAN:

vanDat           <- data.frame(zDat$SalVentAttn, zDat$Estro)
colnames(vanDat) <- c('SalVenAttn', 'Estradiol')

vanVAR <- VAR(vanDat, p = 2)
summary(vanVAR)

vanPerm <- permTestVAR(vanDat, vanVAR, 873)

# Run VAR for SomMot:

smDat           <- data.frame(zDat$SomMot, zDat$Estro)
colnames(smDat) <- c('SomMot', 'Estradiol')

smVAR <- VAR(smDat, p = 2)
summary(smVAR)

smPerm <- permTestVAR(smDat, smVAR, 721)

# Run VAR for TempPar:

tempParDat           <- data.frame(zDat$TempPar, zDat$Estro)
colnames(tempParDat) <- c('TempPar', 'Estradiol')

tempParVAR <- VAR(tempParDat, p = 2)
summary(tempParVAR)

tpPerm <- permTestVAR(tempParDat, tempParVAR, 634)

# Run VAR for Vis:

visDat           <- data.frame(zDat$Vis, zDat$Estro)
colnames(visDat) <- c('Vis', 'Estradiol')

visVAR <- VAR(visDat, p = 2)
summary(visVAR)

visPerm <- permTestVAR(visDat, visVAR, 145)

# Run VAR for Subcort:

scDat           <- data.frame(zDat$Subcort, zDat$Estro)
colnames(scDat) <- c('Subcort', 'Estradiol')

scVAR <- VAR(scDat, p = 2)
summary(scVAR)

scPerm <- permTestVAR(scDat, scVAR, 987)

# Fit models using between-network participation.
#--------------------------------------------------------------------------

# Load in data.

dat <- read.csv('./data.participation.csv')

# Z-transform.

zDat           <- data.frame(scale(dat[,2:27]))
colnames(zDat) <- colnames(dat)[2:27]

# Run VAR for FCN:

fcnDat           <- data.frame(zDat$Control, zDat$Estro)
colnames(fcnDat) <- c('FCN', 'Estradiol')

fcnVAR <- VAR(fcnDat, p = 2)
summary(fcnVAR)

fcnPerm  <- permTestVAR(fcnDat, fcnVAR, 139)

# Run VAR for DMN:

dmnDat           <- data.frame(zDat$DMN, zDat$Estro)
colnames(dmnDat) <- c('DMN', 'Estradiol')

dmnVAR <- VAR(dmnDat, p = 2)
summary(dmnVAR)

dmnPerm <- permTestVAR(dmnDat, dmnVAR, 164)

# Run VAR for DAN:

danDat           <- data.frame(zDat$DorsAttn, zDat$Estro)
colnames(danDat) <- c('DorsAttn', 'Estradiol')

danVAR <- VAR(danDat, p = 2)
summary(danVAR)

danPerm <- permTestVAR(danDat, danVAR, 624)

# Run VAR for Limbic:

limbicDat           <- data.frame(zDat$Limbic, zDat$Estro)
colnames(limbicDat) <- c('Limbic', 'Estradiol')

limbicVAR <- VAR(limbicDat, p = 2)
summary(limbicVAR)

limbicPerm <- permTestVAR(limbicDat, limbicVAR, 852)

# Run VAR for VAN:

vanDat           <- data.frame(zDat$SalVentAttn, zDat$Estro)
colnames(vanDat) <- c('SalVenAttn', 'Estradiol')

vanVAR <- VAR(vanDat, p = 2)
summary(vanVAR)

vanPerm <- permTestVAR(vanDat, vanVAR, 913)

# Run VAR for SomMot:

smDat           <- data.frame(zDat$SomMot, zDat$Estro)
colnames(smDat) <- c('SomMot', 'Estradiol')

smVAR <- VAR(smDat, p = 2)
summary(smVAR)

smPerm <- permTestVAR(smDat, smVAR, 379)

# Run VAR for TempPar:

tempParDat           <- data.frame(zDat$TempPar, zDat$Estro)
colnames(tempParDat) <- c('TempPar', 'Estradiol')

tempParVAR <- VAR(tempParDat, p = 2)
summary(tempParVAR)

tpPerm <- permTestVAR(tempParDat, tempParVAR, 712)

# Run VAR for Vis:

visDat           <- data.frame(zDat$Vis, zDat$Estro)
colnames(visDat) <- c('Vis', 'Estradiol')

visVAR <- VAR(visDat, p = 2)
summary(visVAR)

visPerm <- permTestVAR(visDat, visVAR, 519)

# Run VAR for Subcort:

scDat           <- data.frame(zDat$Subcort, zDat$Estro)
colnames(scDat) <- c('Subcort', 'Estradiol')

scVAR <- VAR(scDat, p = 2)
summary(scVAR)

scPerm <- permTestVAR(scDat, scVAR, 806)
