


sampPreg <- function(dspDat, betaDays, betaCovs, xi=NULL, phi, fwLen, verbose=FALSE) {
  idVec <- unique(dspDat$subjId)
  n <- length(idVec)
  xiVec <- rgamma(n=n, shape=phi, rate=phi)
  beta <- c(betaDays, unlist(betaCovs))
  sexInd <- (dspDat$intercourse == "yes") + 0
  
  # Create U matrix ------------------------------------------------------------
  theModelFormula <- formula( paste("~ -1 + cycleDay + ", paste(names(betaCovs), collapse=" + ")) )
  uMat <- model.matrix(theModelFormula, data=dspDat)
  
  
  # Obtain indices for each subject / cycle ------------------------------------
  idIdx <- lapply(X=idVec, FUN=function(x) which(dspDat$subjId==x))
    # Row indices corresponding to each cycle #
  cycIdx <- lapply(X=idIdx, FUN=function(x) lapply(seq(1,length(x),fwLen), 
                                                     function(j) x[j:(j + fwLen - 1)]) )

 
  # Calculate probability of pregnancy for cycle i,j,k -------------------------
  getPi <- function(thisXi, thisU, thisSex) {
    thisLam <- exp( -1 * thisSex * thisXi * exp( thisU %*% beta ) )
    return ( 1 - prod(thisLam) )
  }
  
   
  # Sample pregnancy for each cycle --------------------------------------------
  #
  # Store as 0 for no, 1 for yes, and 2 for *cycle not needed due to a pregnancy in a
  # previous cycle*
  
  pregVec <- integer( length=(nrow(dspDat) / fwLen) )
  niVec <- sapply(X=cycIdx, FUN=length)
  piVec <- numeric( length=(nrow(dspDat) / fwLen) )
  ctr <- 1
 
  for (i in 1:n) {
    subjPregBool <- FALSE
    
    for (j in 1:niVec[i]) {
      
      if (!subjPregBool) {
        thisIdx <- cycIdx[[ i ]][[ j ]]
        thisU <- uMat[thisIdx, ]
        thisSex <- sexInd[thisIdx]
        
        piVec[ctr] <- getPi( thisXi=xiVec[i], thisU=thisU, thisSex=thisSex )
        pregVec[ctr] <- rbinom(n=1, size=1, prob=piVec[ctr])
        
        if (pregVec[ctr] == 1)
          subjPregBool <- TRUE
      }
      else
        pregVec[ctr] <- 2
      
      ctr <- ctr + 1
    }
  }
  
  if (!verbose)
    return (pregVec)
  else
    return ( list(pregVec, piVec) )
}




# Remove cycles not needed due to successful preg ==============================

rmSuperfluous <- function(baseline, cycle, daily, pregVec, wVec=NULL, xiVec=NULL) {
  fwLen <- nrow(daily) / nrow(cycle)  
  keepBoolCyc <- (pregVec != 2)
  keepBoolDay <- rep(x=keepBoolCyc, each=fwLen)
  
  baseline <- data.frame( subjId = unique(cycle$subjId),
                          baseline )
  cycle <- cycle[keepBoolCyc, ]
  daily <- daily[keepBoolDay, ]
  daily <- data.frame( subjId = rep(cycle$subjId, each=fwLen),
                       cycle = rep(cycle$cycle, each=fwLen),
                       daily )
  pregVec <- pregVec[keepBoolCyc]
  dspDat <- list( baseline=baseline,
                  cycle = data.frame(cycle, pregInd=pregVec),
                  daily = daily)
  if (!is.null(wVec))
    dspDat$W <- wVec[keepBoolDay]
  if (!is.null(xiVec))
    dspDat$xi <- xiVec

  return (dspDat)
}











