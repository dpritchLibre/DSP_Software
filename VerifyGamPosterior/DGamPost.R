
# Density fcn for gamma posterior ----------------------------------------------

dgammaPost <- function(gamVal, gamCoef, gamLoc, W, X, U, xiDay, hypGam) {
  list2env(hypGam, envir=environment())
  betaLvOut <- log(gamCoef[-gamLoc])
  
  aTilde <- getaTilde(W, U, gamLoc, ah)
  bTilde <- getbTilde(gamLoc, X, U, betaLvOut, xiDay, bh)
  pTilde <- getpTilde(ph, ah, bh, bndL, bndU, aTilde, bTilde)

  ptMassBool <- ( gamVal == 1 )
  denVal <- numeric(length=length(gamVal))
  denVal[ptMassBool] <- pTilde
  denVal[!ptMassBool] <- (1 - pTilde) * dTrGamma(x=gamVal[!ptMassBool], a=aTilde, 
                                                 b=bTilde, bndL=bndL, bndU=bndU)
  return (denVal)
}




# Density fcn for gamma posterior Dunson ver -----------------------------------

dgammaPostDuns <- function(gamVal, gamCoef, gamLoc, W, X, U, xiDay, hypGam) {
  list2env(hypGam, envir=environment())
  betaLvOut <- log(gamCoef[-gamLoc])
  
  aTilde <- ah + drop( crossprod(W, U[, gamLoc]) )
  bTilde <- getbTildeDuns(gamLoc, X, U, betaLvOut, xiDay, bh)
  pTilde <- getpTilde(ph, ah, bh, bndL, bndU, aTilde, bTilde)
  
  ptMassBool <- ( gamVal == 1 )
  denVal <- numeric(length=length(gamVal))
  denVal[ptMassBool] <- pTilde
  denVal[!ptMassBool] <- (1 - pTilde) * dTrGamma(x=gamVal[!ptMassBool], a=aTilde, 
                                                 b=bTilde, bndL=bndL, bndU=bndU)
  return (denVal)
}




# Distribution fcn for gamma posterior -----------------------------------------
#
# TODO: Should try to calculate closed-form for this and compare speed

pgammaPost <- function(gamVal, gamCoef, gamLoc, W, X, U, xiDay, p, a, b, bndL, bndU) {
  
  integrate(f=dgammaPost, lower=bndL, upper=gamVal, gamCoef=gamCoef, gamLoc=gamLoc,
            W=W, X=X, U=U, p=p, a=a, b=b, bndL=bndL, bndU=bndU)
}




# Calc the density function of a truncated gamma -------------------------------

dTrGamma <- function(x, a, b, bndL, bndU) {
  dgamma(x=x, shape=a, rate=b) / getTrGamNorm(a, b, bndL, bndU)
}




# Calculate the value of a_h tilde ---------------------------------------------

getaTilde <- function(W, U, gamLoc, ah) {
  ah + drop( crossprod(W, U[, gamLoc]) )
}



# Calculate the value of b_h tilde ---------------------------------------------

getbTilde <- function(gamLoc, X, U, betaLvOut, xiDay, b) {
  
  subsetIdx <- (X == 1) & drop(U[, gamLoc] == 1)
  xiSub <- xiDay[subsetIdx]
  Usub <- U[subsetIdx, -gamLoc]
  logProdTerms <- log(xiSub) + drop(Usub %*% betaLvOut)
  return ( b + sum(exp(logProdTerms)) )
}




# Calculate the value of b_h tilde Dunson ver ----------------------------------
getbTildeDuns <- function(gamLoc, X, U, betaLvOut, xiDay, b) {
  
  subsetIdx <- drop(U[, gamLoc] == 1)
  xiSub <- xiDay[subsetIdx]
  Usub <- U[subsetIdx, -gamLoc]
  logProdTerms <- log(xiSub) + drop(Usub %*% betaLvOut)
  return ( b + sum(exp(logProdTerms)) )
}




# Calculate p_h tilde ----------------------------------------------------------
#
# Defined as p / (p + d2) where d2 is given by:
#
#
#    (1 - p) * C(a,b) * int G(gam; aTilde,bTilde) dgam * exp(bTilde - b)
#    -------------------------------------------------------------------
#                C(aTilde,bTilde) * int G(gam; a,b) dgam
#
#
# Calculations are performed on the log scale

getpTilde <- function(p, a, b, bndL, bndU, aTilde, bTilde) {

  d2numer <- ( log(1 - p) + getLogGamConst(a, b)
               + log( getTrGamNorm(aTilde, bTilde, bndL, bndU) )
               + bTilde - b )
  d2denom <- getLogGamConst(aTilde, bTilde) + log( getTrGamNorm(a, b, bndL, bndU) )
  d2 <- exp( d2numer - d2denom )
  pTilde <- p / (p + d2)
  
  return (pTilde)
}



# Calculate the constant part of a gamma density -------------------------------

getLogGamConst <- function(a, b) {
  a * log(b) - lgamma(a)
}




# Calculate the norming constant for truncated gam -----------------------------

getTrGamNorm <- function(a, b, bndL, bndU) {
  pgamma(q=bndU, shape=a, rate=b) - pgamma(q=bndL, shape=a, rate=b)
}

