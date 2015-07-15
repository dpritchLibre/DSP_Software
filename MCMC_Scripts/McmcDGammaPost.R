
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
# Calculations are performed on the log scale

getpTilde <- function(p, a, b, bndL, bndU, aTilde, bTilde) {

  d2numer <- ( log(1 - p) + getLogGamConst(a, b)
               + log( getNormingConst(aTilde, bTilde, bndL, bndU) )
               + bTilde - b )
  d2denom <- getLogGamConst(aTilde, bTilde) + log( getNormingConst(a, b, bndL, bndU) )
  d2 <- exp( d2numer - d2denom )
  pTilde <- p / (p + d2)
  
  return (pTilde)
}



# Calculate the constant part of a gamma density -------------------------------

getLogGamConst <- function(a, b) {
  a * log(b) - lgamma(a)
}




# Calculate the norming constant for truncated gam -----------------------------

getNormingConst <- function(a, b, bndL, bndU) {
  pgamma(q=bndU, shape=a, rate=b) - pgamma(q=bndL, shape=a, rate=b)
}









  


