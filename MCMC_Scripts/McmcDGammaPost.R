
# Density fcn for gamma posterior ----------------------------------------------

dgammaPost <- function(gamVal, gamCoef, gamLoc, W, X, U, p, a, b, bndL, bndU) {
  #cat("gamVal: ", gamVal, "\n")
  betaLvOut <- log(gamCoef)
  aTilde <- p + sum(W * U[, gamLoc] )
  
  subsetIdx <- (X == 1) & drop(U[, gamLoc] == 1)
  xiSub <- xiDay[subsetIdx]
  Usub <- U[subsetIdx, -gamLoc]
  logProdTerms <- log(xiSub[x]) + drop(Usub %*% betaLvOut)
  bTilde <- b + sum(exp(logProdTerms))
  
  d2numer <- ( (1 - p) * getGammaConst(a, b) * getNormingConst(aTilde, bTilde, bndL, bndU) 
               * exp(bTilde - b) )
  d2denom <- getGammaConst(aTilde, bTilde) * getNormingConst(a, b, bndL, bndU)
  d2 <- d2numer / d2denom
  pTilde <- p / (p + d2)

  ptMassBool <- ( gamVal == 1 )
  denVal <- numeric(length=length(gamVal))
  denVal[ptMassBool] <- pTilde
  denVal[!ptMassBool] <- (1 - pTilde) * dTrGamma(x=gamVal[!ptMassBool], a=aTilde, 
                                                 b=bTilde, bndL=bndL, bndU=bndU)
  
  return (denVal)
}




# Distribution fcn for gamma posterior -----------------------------------------
#
# Should try to calculate closed-form for this

pgammaPost <- function(gamVal, gamCoef, gamLoc, W, X, U, p, a, b, bndL, bndU) {
  
  integrate(f=dgammaPost, lower=bndL, upper=gamVal, gamCoef=gamCoef, gamLoc=gamLoc,
            W=W, X=X, U=U, p=p, a=a, b=b, bndL=bndL, bndU=bndU)
}



# Calc the density function of a truncated gamma -------------------------------

dTrGamma <- function(x, a, b, bndL, bndU) {
  
  return (dgamma(x=x, shape=a, rate=b) / getNormingConst(a, b, bndL, bndU))
}



# Calculate the constant part of a gamma density -------------------------------
getGammaConst <- function(a, b) {
  
  return ( b^a / gamma(a) )  
}




# Calculate the norming constant for truncated gam -----------------------------
getNormingConst <- function(a, b, bndL, bndU) {
  
  return ( pgamma(q=bndU, shape=a, rate=b) - pgamma(q=bndL, shape=a, rate=b) )
}




  


