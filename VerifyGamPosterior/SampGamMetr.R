
sampGamMetrop <- function(W, X, U, gamCoef, xiDay, gamLoc, gamInit, 
                          hypGam, numSamp, pJ, delta, trackProg=TRUE) {
  currTheta <- gamInit
  betaLvOut <- log(gamCoef[-gamLoc])
  gamSampVec <- numeric(numSamp)
  trackVals <- sapply(seq(0.1, 0.9, by=0.1), function(x) max(which( 1:numSamp / numSamp <= x)))

  for (i in 1:numSamp) {
    
    # Sample M
    currM <- sampM(W, X, U, gamCoef, xiDay, gamLoc, currTheta, hypGam)
    
    # Sample theta_h0
    if (currM) {
      currTheta <- rgamma(1, shape=hypGam$ah, rate=hypGam$bh)
      gamSampVec[i] <- 1
    }
    else {
      aTilde <- getaTilde(W, U, gamLoc, hypGam$ah)
      bTilde <- getbTilde(gamLoc, X, U, betaLvOut, xiDay, hypGam$bh)
      currTheta <- rgamma(1, shape=aTilde, rate=bTilde)
      gamSampVec[i] <- currTheta
    }
    
    gamCoef[gamLoc] <- gamSampVec[i]
    if (trackProg && (i %in% trackVals))
      cat(round(100 * i / numSamp), "%..  ", sep="")
  }
  if (trackProg) cat("\n")
  
  return (gamSampVec)
}




# Calculate posterior prob that M == 1 -----------------------------------------

sampM <- function(W, X, U, gamCoef, xiDay, gamLoc, currTheta, hypGam) {
  gamIsOneBeta <- replace(log(gamCoef), list=gamLoc, values=0)
  gamIsNotOneBeta <- replace(log(gamCoef), list=gamLoc, values=log(currTheta))
  
  termGamIsOne <- dW(W, X, U, gamIsOneBeta, xiDay, FALSE) * hypGam$ph
  termGamIsNotOne <- dW(W, X, U, gamIsNotOneBeta, xiDay, FALSE) * (1 - hypGam$ph)
  
  mProb <- termGamIsOne / (termGamIsOne + termGamIsNotOne)
  
  return ( as.logical( rbinom(1, size=1, prob=mProb) ) )
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



