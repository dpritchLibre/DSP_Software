
# Sample xi from conditional dist ----------------------------------------------

sampXi <- function(W, uProdBeta, phi, idIdx) {
  
  a <- phi + sapply(X=idIdx, FUN=function(x) sum(W[x]))
  b <- sapply(X=idIdx, FUN=function(x) sum(exp(uProdBeta[x])))
  
  return ( rgamma(n=length(idIdx), shape=a, rate=b) )
}




# Sample phi from conditional dist ---------------------------------------------

sampPhiProp <- function(phi, phiPropDist="uniform", delta) {
  
  if (phiPropDist == "normal")
    phiProp <- abs( rnorm(n=1, mean=phi, sd=delta) )
  else
    phiProp <- abs( runif(n=1, min=(phi - delta), max=(phi + delta)) )
  
  return (phiProp)
}




# Calculate log(r) -------------------------------------------------------------
#
# r is given by
#
#           ( prod_i pi(xi_i | phiProposal) ) * pi(phiProposal)
#           ---------------------------------------------------
#            ( prod_i pi(xi_i | phiCurrent) ) * pi(phiCurrent)

getPhiLogR <- function(xi, phiCurr, phiProp, hypPhi) {
  
  logNumerTerm1 <- sum( dgamma(x=xi, shape=phiProp, rate=phiProp, log=TRUE) )
  logNumerTerm2 <- log( dgamma(x=phiProp, shape=hypPhi$c1, rate=hypPhi$c2) )
  
  logDenomTerm1 <- sum( dgamma(x=xi, shape=phiCurr, rate=phiCurr, log=TRUE ) )
  logDenomTerm2 <- log( dgamma(x=phiCurr, shape=hypPhi$c1, rate=hypPhi$c2) )
  
  return( logNumerTerm1 + logNumerTerm2 - logDenomTerm1 - logDenomTerm2 )
}




