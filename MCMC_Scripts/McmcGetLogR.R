
# Calculate the acceptance ratio r ---------------------------------------------

getLogR <- function(Y, U, xi, betaCoef, currVal, proposeVal, changeLoc, hyperparGam) {
  
  propW <- calcLik_wCondGam(Y=Y, U=U, xi=xi, betaCoef=betaCoef, 
                            updateVal=proposeVal, changeLoc=changeLoc)
  propPhi <- calcLik_gammaH(updateVal=proposeVal, hyperparGam=hyperparGam)
  currW <- calcLik_wCondGam(Y=Y, U=U, xi=xi, betaCoef=betaCoef, 
                            updateVal=currVal, changeLoc=changeLoc)
  currPhi <- calcLik_gammaH(updateVal=currVal, hyperparGam=hyperparGam)
  
  return (propW + propPhi - currW - currPhi)
}




# Calculate log-lik value for p(W | gamma, ...) --------------------------------

calcLik_wCondGam <- function(Y, U, xi, betaCoef, updateVal, changeLoc) {
  betaCoef[changeLoc] <- updateVal
  
  # Days for which conception occured during cycle #
  gotPregBool <- ( rep(x=Y, each=fwLen) == 1 )
  meanVec <- xi[gotPregBool] * drop( exp(U[gotPregBool, ] %*% betaCoef) )
  
  logProd <- sum(log( dpois(x=W[gotPregBool], lambda=meanVec) ))
  return (logProd)
}




# Calculate log-lik value for p(gamma_h) ---------------------------------------

calcLik_gammaH <- function(updateVal, hyperparGam) {
  updateGam <- exp(updateVal)
  p <- hyperparGam$ptMassProp
  a <- hyperparGam$shape
  b <- hyperparGam$rate
  bndL <- hyperparGam$boundL
  bndU <- hyperparGam$boundU
  
  if (updateGam == 1)
    return ( log(p) )
  else {
    normingVal <- pgamma(q=bndU, shape=a, rate=b) - pgamma(q=bndL, shape=a, rate=b)
    dTrunGamma <- dgamma(x=updateGam, shape=a, rate=b) / normingVal
    return ( log( (1 - p) * dTrunGamma ) ) 
  }
}