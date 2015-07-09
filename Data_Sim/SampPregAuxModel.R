
sampPregAuxModel <- function(dspDat, betaDays, betaCovs, xi=NULL, phi, fwLen, verbose=FALSE) {
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
  
  
  # Calculate lambda parameter for each day in cycle day i,j -------------------
  getLam <- function(thisXi, thisU, thisSex) {
    thisLam <- thisSex * thisXi * exp( thisU %*% beta )
    return (thisLam)
  }
  
  
  # Sample pregnancy for each cycle --------------------------------------------
  #
  # Store as 0 for no, 1 for yes, and 2 for *cycle not needed due to a pregnancy in a
  # previous cycle*
  
  pregVec <- integer(length=(nrow(dspDat) / fwLen))
  niVec <- sapply(X=cycIdx, FUN=length)
  lamVec <- numeric(length=nrow(dspDat))
  wVec <- integer(length=nrow(dspDat))
  ctr <- 1
  
  for (i in 1:n) {
    subjPregBool <- FALSE
    
    for (j in 1:niVec[i]) {
      
      if (!subjPregBool) {
        thisIdx <- cycIdx[[ i ]][[ j ]]
        thisU <- uMat[thisIdx, ]
        thisSex <- sexInd[thisIdx]
        
        lamVec[thisIdx] <- getLam( thisXi=xiVec[i], thisU=thisU, thisSex=thisSex )
        wVec[thisIdx] <- rpois(n=fwLen, lambda=lamVec[thisIdx])
        if (sum(wVec[thisIdx]) > 0) {
          pregVec[ctr] <- 1
          subjPregBool <- TRUE
        }
      }
      else
        pregVec[ctr] <- 2
      
      ctr <- ctr + 1
    }
  }
  
  if (!verbose)
    return (pregVec)
  else
    return ( list(pregVecPre=pregVec, lamVecPre=lamVec, wVecPre=wVec, xiVecPre=xiVec) )
}