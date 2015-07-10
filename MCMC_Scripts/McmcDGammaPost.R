
# Density fcn for gamma posterior ----------------------------------------------

dgammaPost <- function(gamVal, gamCoef, gamLoc, W, X, U, xi, p, a, b, bndL, bndU) {
  betaLvOut <- log(gamCoef[-gamLoc])
  
  aTilde <- a + sum(W * U[, gamLoc] )
  bTilde <- getbTilde(gamLoc, X, U, xi, b, betaLvOut)
  pTilde <- getpTilde(p, a, b, bndL, bndU, aTilde, bTilde)

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

pgammaPost <- function(gamVal, gamCoef, gamLoc, W, X, U, xi, p, a, b, bndL, bndU) {
  
  integrate(f=dgammaPost, lower=bndL, upper=gamVal, gamCoef=gamCoef, gamLoc=gamLoc,
            W=W, X=X, U=U, p=p, a=a, b=b, bndL=bndL, bndU=bndU)
}




# Calc the density function of a truncated gamma -------------------------------

dTrGamma <- function(x, a, b, bndL, bndU) {
  
  return (dgamma(x=x, shape=a, rate=b) / getNormingConst(a, b, bndL, bndU))
}



# Calculate the value of a_h tilde  --------------------------------------------

getbTilde <- function(gamLoc, X, U, xi, b, betaLvOut) {
  
  subsetIdx <- (X == 1) & drop(U[, gamLoc] == 1)
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

getpTilde <- function(p, a, b, bndL, bndU, aTilde, bTilde) {

  d2numer <- ( (1 - p) * getGammaConst(a, b) * getNormingConst(aTilde, bTilde, bndL, bndU) 
               * exp(bTilde - b) )
  d2denom <- getGammaConst(aTilde, bTilde) * getNormingConst(a, b, bndL, bndU)
  d2 <- d2numer / d2denom
  pTilde <- p / (p + d2)
}



# Calculate the constant part of a gamma density -------------------------------

getGammaConst <- function(a, b) {

  suppressWarnings( try(return( b^a / gamma(b))) )
  return (0)  
}




# Calculate the norming constant for truncated gam -----------------------------
#
# TODO: how does this behave for large aTilde, bTilde when bounds not 0, 1?

getNormingConst <- function(a, b, bndL, bndU) {
  
  return ( pgamma(q=bndU, shape=a, rate=b) - pgamma(q=bndL, shape=a, rate=b) )
}









  


