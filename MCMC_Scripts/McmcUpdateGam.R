
# Sample gamma_h from full conditional posterior ------------------------------

sampGammaH <- function(W, uProdBetaNoH, xiDay, hypGamH, uIsOneBool, pregUIsOneBool, hIsTrunBool) {
  list2env(hypGamH, envir=environment())
  
  aTilde <- a + sum(W[pregUIsOneBool])
  bTilde <- getbTilde(uProdBetaNoH, xiDay, b, uIsOneBool)
  pTilde <- getpTilde(p, a, b, bndL, bndU, aTilde, bTilde)
  
  if (sampBool(pTilde))
    return (1)
  else
    return ( sampGammaTr(aTilde, bTilde, bndL, bndU, hIsTrunBool) )
}




# Sample a scalar from a truncated gamma dist ----------------------------------

sampGammaTr <- function(aTilde, bTilde, bndL, bndU, hIsTrunBool) {
  
  if (!hIsTrunBool)
    gamSampVal <- rgamma(n=1, shape=aTilde, rate=bTilde)
  else {
    unifBndL <- pgamma(bndL, shape=aTilde, rate=bTilde)
    unifBndU <- pgamma(bndU, shape=aTilde, rate=bTilde)
    unifSampVal <- runif(n=1, min=unifBndL, max=unifBndU)
    gamSampVal <- qgamma(unifSampVal, shape=aTilde, rate=bTilde)
  }

  return (gamSampVal)
}




# Calculate the value of b_h tilde  --------------------------------------------
#
# Defined as 
#
#    bh + sum_{i,j,k: X_ijk == 1, u_ijkh == 1} ( xi_i * prod_{l ne h} gamma_l^(u_ijkl) )
#
# Calculates the terms inside the summation on the log scale

getbTilde <- function(uProdBetaNoH, xiDay, b, uIsOneBool) {

  logProdVec <- log( xiDay[uIsOneBool] ) + uProdBetaNoH[uIsOneBool]
  return ( b + sum(exp(logProdVec)) )
}




# Calculate p_h tilde ----------------------------------------------------------
#
# Defined as p / (p + d2) where d2 is given by:
#
#
#    (1 - p) * C(a,b) * int{ G(gam; aTilde,bTilde) }dgam * exp(bTilde - b)
#    -------------------------------------------------------------------
#                C(aTilde,bTilde) * int{ G(gam; a,b) }dgam
#
#
# Calculations are performed on the log scale

getpTilde <- function(p, a, b, bndL, bndU, aTilde, bTilde) {
  
  d2numer <- ( log(1 - p) 
               + getLogGamConst(a, b)
               + getLogTrGamNorm(aTilde, bTilde, bndL, bndU)
               + bTilde
               - b )
  d2denom <- getLogGamConst(aTilde, bTilde) + getLogTrGamNorm(a, b, bndL, bndU)
  d2 <- exp( d2numer - d2denom )
  pTilde <- p / (p + d2)
  
  return (pTilde)
}




# Calculate the constant term from a gamma density -------------------------------

getLogGamConst <- function(a, b) {
  a * log(b) - lgamma(a)
}




# Calculate the norming constant for truncated gam -----------------------------

getLogTrGamNorm <- function(a, b, bndL, bndU) {
  log( pgamma(q=bndU, shape=a, rate=b) - pgamma(q=bndL, shape=a, rate=b) )
}




# Sample Boolean with probability 'prob' for TRUE ------------------------------

sampBool <- function(prob) {
  sample(c(TRUE, FALSE), size=1, prob=c(prob, 1 - prob))
}