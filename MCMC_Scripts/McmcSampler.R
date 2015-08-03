
# Read in data -----------------------------------------------------------------
{
  if (Sys.info()[8] == "dpritch") {
    setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
    outPath <- "/home/dpritch/Documents/Projects/Dunson Day Specific/Software/Output/"
  }
  else {
    setwd("/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/")
    outPath <- "/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/Output/"
  }
}

source("MCMC_Scripts/McmcHelperW.R")
source("MCMC_Scripts/McmcHelperGam.R")
source("MCMC_Scripts/McmcHelperXi.R")
source("MCMC_Scripts/McmcHelperPhi.R")
source("Format_Data/FormatDspDat.R")
load("Data/RawData.RData")





# Prepare data for use in MCMC sampler -----------------------------------------

idName <- "subjId"
cycName <- "cycle"
pregName <- "pregInd"
sexName <- "intercourse"
fwInd <- "fwInd"
varInclNames <- list( baseline = c("age","bmi","gravid"),
                      cycle = "cycleLen",
                      daily = "lube" )
fwLen <- 5

# Delete some arbitrary data to make sure that program can put the data together properly
myDspDat <- dspDat(baseline[-c(1,2), ], cycle[-(20:25), ], daily[-tail(1:nrow(daily), n=5), ], 
                   idName, cycName, pregName, sexName, fwInd, varInclNames, fwLen)
summary(myDspDat)

# Place the sampler objects into global environment for running script
list2env(myDspDat$samplerObj, envir=environment())




# Set hyperparams --------------------------------------------------------------

###Initial values
phi <- 1
gamCoef <- rep(1, q)
uProdBeta <- drop( U %*% log(gamCoef) )
xi <- rep(1, n)
xiDay <- xi[idDayExpan]

###Hyperparameters
hypGam <- list(ah=1, bh=1, ph=0.5, bndL=0, bndU=5)
hypPhi <- list(c1=1, c2=1)

###Metropolis objects
delta <- 0.2
numAcceptPhi <- 0

###MCMC Objects
numSamp <- 5000
write(1:q, file=paste(outPath,"GAMMA.csv",sep=""), sep=",", ncolumns=q)
write(1:n, file=paste(outPath,"XI.csv",sep=""), sep=",", ncolumns=n)
write(1, file=paste(outPath,"PHI.csv",sep=""), sep=",", ncolumns=1)

###Time MCMC Sampler
#begin <- Sys.time()














# MCMC Sampler =================================================================

for (s in 1:numSamp) {
  
  # Sample latent variable W
  W <- sampW(uProdBeta, xiDay, pregDayBool, pregCycIdx)

  # Sample regression coefficients gamma
  for (h in 1:q) {
    uProdBetaNoH <- uProdBeta - drop(U[, h] * log(gamCoef[h]))
    gamCoef[h] <- sampGammaH(W, uProdBetaNoH, xiDay, hypGam, uBool[[h]], 
                             pregUBool[[h]], TRUE)
    uProdBeta <- uProdBetaNoH + drop(U[, h] * log(gamCoef[h]))
  }

  # Sample woman-specific fecundability multiplier xi
  xi <- sampXi(W, uProdBeta, phi, idIdx, idPregIdx, pregCycIdx, n)
  xiDay <- xi[idDayExpan]
  
  # Metroplois step for phi, the variance parameter for xi
  phiProp <- sampPhiProp(phi, delta)
  phiLogR <- getPhiLogR(xi=xi, phiCurr=phi, phiProp=phiProp, hypPhi=hypPhi)
  if (log(runif(1)) < phiLogR) {
    phi <- phiProp
    numAcceptPhi <- numAcceptPhi + 1
  }
  
  # Write samples to output
  write(phi, file=paste0(outPath, "PHI.csv"), sep=",", ncolumns=1, append=TRUE)
  write(xi, file=paste0(outPath, "XI.csv"), sep=",", ncolumns=n, append=TRUE)
  write(gamCoef, file=paste0(outPath, "GAMMA.csv"), sep=",", ncolumns=q, append=TRUE)
  
  if ((s %% 1000) == 0)
    print(s)
  
#   # Verbose (we can personalize the verbose that we output...)
#   if (s!=nsims) print(paste("Completed Percentage: ",round((s/nsims)*100,digits=0),"%",sep=""))
#   print(paste("Gamma: ",round(gamma[1],digits=3),", ",round(gamma[2],digits=3),sep=""))
#   print(paste("Phi: ",round(phi,digits=3),sep=""))
#   print(paste("Deviance: ",round(dev,digits=3),sep=""))
#   print(paste("Acceptance Rate: Phi (",100*round(mean(acceptance_phi)/s,digits=2),"%)",sep=""))
#   print("######################################################################################################")
#   ###Time MCMC Sampler
#   if (s==nsims) {
#     after<-Sys.time()
#     time<-after-begin
#     print(paste("Run Time: ",round(time,digits=2)," ",attr(time,"unit"),sep=""))
#     print("######################################################################################################")
#   }
  
  ###End MCMC Sampler	
}





# Summary stats of the output --------------------------------------------------

# Gamma coefs
gamTab <- apply(read.csv(file=paste0(outPath, "GAMMA.csv")), MARGIN=2, 
                FUN=quantile, probs=c(0.025, 0.500, 0.975))
trueVals <- c(0.14, 0.08, 0.34, 0.31, 0.08, exp(c(-0.08, -0.43, -1.03, -0.22, -0.47, 1.21)), rep(1,3))
rbind(gamTab, trueVals)

