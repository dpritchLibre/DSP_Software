
sampGamMetrop <- function(W, X, U, gamCoef, xiDay, gamLoc, gamInit, 
                          hypGam, B, pJ, delta, verbose=FALSE, trackProg=TRUE) {
  betaCoef <- log(gamCoef)
  currGam <- gamInit
  gamSampVec <- numeric(length=B)
  if (verbose)
    logRVec <- propValVec <- numeric(length=B)
  gamCurrCoef <- replace(x=gamCoef, list=gamLoc, values=gamInit)
  trackVals <- sapply(X=seq(0.1, 0.9, by=0.1), FUN=function(x) max(which( 1:B / B <= x)))
 
  numAccept <- 0
  sampBool <- function(p) sample(c(TRUE, FALSE), size=1, prob=c(p, 1-p))
  
  
  for (i in 1:B) {
    
    if (sampBool(pJ))
      propVal <- 1
    else 
      propVal <- abs( runif(n=1, min=(currGam - delta), max=(currGam + delta)) )

    logR <- getLogR(W, X, U, gamCurrCoef, xiDay, propVal, gamLoc, hypGam)
    if (log(runif(1)) < logR) {
      gamCurrCoef[gamLoc] <- currGam <- propVal
      numAccept <- numAccept + 1
    }

    gamSampVec[i] <- currGam
    if (verbose) {
      logRVec[i] <- logR
      propValVec[i] <- propVal
    }
    
    if (trackProg && (i %in% trackVals))
      cat(round(100 * i / B), "%..  ", sep="")
  }
  if (trackProg)
    cat("\n")
  
  if (!verbose) 
    return (gamSampVec)
  else {
    gamMetropList <- list( gamSamp = gamSampVec,
                           acceptRate = numAccept / B,
                           propVal = propValVec,
                           logR = logRVec )
    return (gamMetropList)
  }
}



# Calculate log gamma acceptance ratio r ---------------------------------------
#
# Acceptance ratio for gamma_h is given by
#
#
#            p(W | gamma*, xi, data) * p(gamma_h*)
#         -------------------------------------------
#         p(W | gamma^(s), xi, data) * p(gamma_h^(s))
# 
#
# where gamma* denotes the gamma vector with the h-th term replaced by the proposal 
# value and similarly for gamma^(s)

getLogR <- function(W, X, U, gamCurrCoef, xiDay, propVal, propLoc, hypGam) {
  betaCoef <- log(gamCurrCoef)
  betaProp <- replace(x=betaCoef, list=propLoc, values=log(propVal))
  wLogLik <- list()
  gammaHLogLik <- list()
  
  wLogLik$prop <- dW(W=W, X=X, U=U, betaCoef=betaProp, xiDay=xiDay)
  gammaHLogLik$prop <- dgammaH(gamVal=propVal, hypGam=hypGam)
  
  wLogLik$curr <- dW(W=W, X=X, U=U, betaCoef=betaCoef, xiDay=xiDay)
  gammaHLogLik$curr <- dgammaH(gamVal=gamCoef[propLoc], hypGam=hypGam)
  
  return (wLogLik$prop + gammaHLogLik$prop - wLogLik$curr - gammaHLogLik$curr)
}




# Calculate log-lik value for p(W | gamma, xi) --------------------------------

dW <- function(W, X, U, betaCoef, xiDay, logLik=TRUE) {

  hadSexBool <- (X == 1)
  meanVec <- xiDay[hadSexBool] * drop( exp( U[hadSexBool, ] %*% betaCoef ) )
  logProd <- sum( dpois(x=W[hadSexBool], lambda=meanVec, log=TRUE) )
  
  if (logLik)
    return (logProd)
  else
    return (exp(logProd))
}




# Calculate log-lik value for p(gamma_h) ---------------------------------------

dgammaH <- function(gamVal, hypGam, logLik=TRUE) {
  list2env(hypGam, envir=environment())
  
  if (gamVal == 1)
    logLikVal <- log(ph)
  else {
    logNrmVal <- log( pgamma(q=bndU, shape=ah, rate=bh) - pgamma(q=bndL, shape=ah, rate=bh) )
    logTrGamma <- dgamma(x=gamVal, shape=ah, rate=bh, log=TRUE) - logNrmVal
    logLikVal <- log(1 - ph) + logTrGamma
  }
  
  if (logLik)
    return (logLikVal)
  else
    return ( exp(logLikVal) )
}



