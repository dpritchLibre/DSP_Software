
dsp <- function(dspDat, numSamp=1e4, hypGam=NULL, hypPhi=NULL, tuningPhi=0.1, tuningGam=NULL,
                trackProg=TRUE, progQuants=seq(0.1, 0.9, 0.1), verbose=FALSE, saveToFile=FALSE, 
                keepOnlyGam=FALSE) {
  # TODO: check if valid input
  source("MCMC_Scripts/McmcHelper.R")
  source("MCMC_Scripts/McmcUpdateW.R")
  source("MCMC_Scripts/McmcUpdateGam.R")
  source("MCMC_Scripts/McmcUpdateXi.R")
  source("MCMC_Scripts/McmcUpdatePhi.R")
  list2env(dspDat$samplerObj, envir=environment())
  
  #outPath <- format_outPath(outPath) # paste0(outPath, "/")
  if (trackProg)
    cat("Progress:  ")
  if (trackProg || verbose) {
    # Sampler iterations for which we print the percentage of progress
    trackVals <- sapply(progQuants, function(x) tail(which(1:numSamp <= x * numSamp), 1))
  }
  
  # Combine default hyperparameters with custom user input hyperparameters
  hypPhi <- getHypPhi(hypPhi)
  hypGam <- getHypGam(varNames, hypGam)
  gamIsTrunBool <- sapply(1:q, function(j) with((a == 0) && (b == Inf), data=hypGam[[j]]))
  
  # Set initial values; uses mean of prior dists for phi and gamma
  phi <- hypPhi$c1 / hypPhi$c2
  gamCoef <- getGamInit(hypGam, gamIsTrunBool)
  uProdBeta <- drop( U %*% log(gamCoef) )
  xi <- rep(1, n)
  xiDay <- xi[idDayExpan]
  
  # Metropolis acceptance rate counters
  numAcceptPhi <- 0
  numAcceptGam <- NULL # TODO: create initialization fcn for continuous gammas
  
  # Inititalize MCMC output files or objects
  if (saveToFile) {
    write(varNames, file=paste0(outPath, "GAMMA.csv"), sep=",", ncolumns=q) # need to update 'getSamplerObj'
    write(subjId, file=paste0(outPath, "XI.csv"), sep=",", ncolumns=n) # need to update 'getSamplerObj'
    write("phi", file=paste0(outPath, "PHI.csv"), sep=",", ncolumns=1)
  }
  else {
    phiOut <- numeric(numSamp)
    xiOut <- data.frame( matrix(nrow=numSamp, ncol=n, dimnames=list(1:numSamp, subjId)) )
    gamOut <- data.frame( matrix(nrow=numSamp, ncol=q, dimnames=list(1:numSamp, varNames)) )
  }

  
  
  
  # MCMC Sampler =================================================================
  
  for (s in 1:numSamp) {
    
    # Sample latent variable W
    W <- sampW(uProdBeta, xiDay, pregDayBool, pregCycIdx)
  
    # Sample regression coefficients gamma
    for (h in 1:q) {
      uProdBetaNoH <- uProdBeta - drop(U[, h] * log(gamCoef[h]))
      gamCoef[h] <- sampGammaH(W, uProdBetaNoH, xiDay, 
                               hypGam[[h]], uBool[[h]], pregUBool[[h]], TRUE)
      uProdBeta <- uProdBetaNoH + drop(U[, h] * log(gamCoef[h]))
    }
  
    # Sample woman-specific fecundability multiplier xi
    xi <- sampXi(W, uProdBeta, phi, idIdx, idPregIdx, pregCycIdx, n)
    xiDay <- xi[idDayExpan]
    
    # Metroplois step for phi, the variance parameter for xi
    phiProp <- sampPhiProp(phi, tuningPhi)
    phiLogR <- getPhiLogR(xi=xi, phiCurr=phi, phiProp=phiProp, hypPhi=hypPhi)
    if (log(runif(1)) < phiLogR) {
      phi <- phiProp
      numAcceptPhi <- numAcceptPhi + 1
    }
    
    # Write samples to output
    if (saveToFile) {
      write(phi, file=paste0(outPath, "PHI.csv"), sep=",", ncolumns=1, append=TRUE)
      write(xi, file=paste0(outPath, "XI.csv"), sep=",", ncolumns=n, append=TRUE)
      write(gamCoef, file=paste0(outPath, "GAMMA.csv"), sep=",", ncolumns=q, append=TRUE)
    }
    else {
      phiOut[s] <- phi
      xiOut[s, ] <- xi
      gamOut[s, ] <- gamCoef
    }
    
    if (trackProg && (s %in% trackVals))
      cat(round(100 * s / numSamp), "%..  ", sep="")
    else if (verbose && (s %in% trackVals)) {
      # verbose output  
    }
    
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

  if (!saveToFile) {
    return ( list( phi = phiOut,
                   xi  = xiOut,
                   gam = gamOut ) )
  }
}

