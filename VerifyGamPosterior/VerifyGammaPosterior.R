
# Read in data -----------------------------------------------------------------
{
  if (Sys.info()[8] == "dpritch")
    setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
  else
    setwd("/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/")
}

load("Data/PracticeDatAuxMod.RData")
source("VerifyGamPosterior/SampGamMetr.R")
source("VerifyGamPosterior/DGamPost.R")
set.seed(0)




# Dataset attributes -----------------------------------------------------------

fwLen <- 5
U <- as.matrix( U[, 1:11] )  # Only keep terms with nonzero beta coefs
Yday <- rep(Y, each=fwLen)




# Set (true) parameter values --------------------------------------------------

betaDays <- log( c(0.14, 0.08, 0.34, 0.31, 0.08) )
betaCovs <- list( age = c(-0.08, -0.43, -1.03),
                  bmi = c(-0.22, -0.47),
                  gravid = 1.21 )
betaCoef <- c(betaDays, unlist(betaCovs))
names(betaCoef)[1:5] <- paste0("day", 1:5)
gamCoef <- exp(betaCoef)

phi <- 1
xiDay <- rep(xi, times=table(id))




# Sample p(gamma_h | ...) via Metropolis algorithm =================================

B <- 1e3        # number of Metropolis iterations
pJ <- 0.1       # prob of selecting gamma_h = 1 in proposal dist
delta <- 0.1    # proposal dist tuning param
hypGam <- list(bndL=0, bndU=Inf, ph=0.5, ah=1, bh=1)


# Sample B gamma values for each gamma_h
gamMetrMat <- matrix(nrow=B, ncol=length(gamCoef))

for (j in 1:length(gamCoef)) {
  cat("Sampling gamma_", j, ":  ", sep="")
  gamInit <- abs( runif(n=1, min=(gamCoef[j] - delta), max=(gamCoef[j] + delta)) )
  gamMetrMat[,j] <- sampGamMetrop(W, X, U, gamCoef, xiDay, j, gamInit, hypGam, B, pJ, delta)
}




# Empirical density from Metropolis samples fcn --------------------------------
#
# 'domain' is a length-2 vector.  Only returns the continuous part of the density

getEmpDen <- function(x, domain, numBins) {
  xCont <- x[abs(x - 1) > 0.0001]
  bins <- seq(from=domain[1], to=domain[2], length.out=(numBins + 1))
  binL <- bins[-(numBins + 1)]
  binU <- bins[-1]
  binCount <- sapply(X=1:numBins, FUN=function(i) sum((binL[i] < xCont) & (xCont <= binU[i])))
  
  wPerBin <- (domain[2] - domain[1]) / numBins
  return (binCount / (length(x) * wPerBin))
}




# Plot Metrop and true density fcn ---------------------------------------------
#
# 'gamMetrSamps' is a vector of Metropolis samples.  'domain' is a length-2 vector

plotGamMetrDen <- function (gamMetrVec, gamLoc, domain, numPts, trueVal, showDuns=FALSE) {
  xVals <- seq(from=domain[1], to=domain[2], length.out=(numPts + 1))[-1]
  
  # Get true posterior density values
  yTrue <- dgammaPost(gamVal=xVals, gamCoef=gamCoef, gamLoc=gamLoc, 
                      W=W, X=X, U=U, xiDay=xiDay, hypGam=hypGam)
  yTrueMass <- dgammaPost(gamVal=1, gamCoef=gamCoef, gamLoc=gamLoc, 
                          W=W, X=X, U=U, xiDay=xiDay, hypGam=hypGam)
  
  # Get Metropolis empirical density values
  yMetr <- getEmpDen(x=gamMetrVec, domain=domain, numBins=numPts)
  yMetrMass <- sum(abs(gamMetrVec - 1) < 0.0001) / length(gamMetrVec)
  
  # Get Dunson posterior density values
  yDuns <- dgammaPostDuns(gamVal=xVals, gamCoef=gamCoef, gamLoc=gamLoc, 
                          W=W, X=X, U=U, xiDay=xiDay, hypGam=hypGam)
  yDunsMass <- dgammaPostDuns(gamVal=1, gamCoef=gamCoef, gamLoc=gamLoc, 
                              W=W, X=X, U=U, xiDay=xiDay, hypGam=hypGam)
  
  # ```````````````` #
  #  Construct plot  #
  # ................ #
  
  par(mar=c(2.5, 2.5, 0.1, 0.1), las=1)
  
  # True density plot
  plot(x=xVals, y=yTrue, ylim=c(0, max(yTrue, yMetr)), ann=FALSE)
  if (yTrueMass > 0.01)
    points(x=1, y=yTrueMass, pch=16)
  
  # Metropolis empirical density plot
  points(x=xVals, y=yMetr, pch=4, col="red", cex=0.75)
  if (yMetrMass > 0.01)
    points(x=1, y=yMetrMass, pch=16, col="red")
  
  # Dunson density plot
  if (showDuns) {
    points(x=xVals, y=yDuns, pch=2, col="green", cex=0.75)
    if (yMetrMass > 0.01)
      points(x=1, y=yDunsMass, pch=16, col="green")
    
    ptMassText <- c(paste0("Point mass true: ", round(yTrueMass, 2)), 
                    paste0("Point mass Metr: ", round(yMetrMass, 2)),
                    paste0("Point mass Duns: ", round(yDunsMass, 2)))
    legend(x="topright", legend=ptMassText, text.col=c("black","red","green"), cex=0.75)
  }
  else {
    ptMassText <- c(paste0("Point mass true: ", round(yTrueMass, 2)), 
                    paste0("Point mass Metr: ", round(yMetrMass, 2)))
    legend(x="topright", legend=ptMassText, text.col=c("black","red"), cex=0.75)
  }
  
  abline(v=trueVal, col="grey50", lwd=2, lty=2)
}



# Convenience wrapper for plotGamMetrDen fcn -----------------------------------

plotGamPost <- function(gamLoc, showDuns=FALSE) {
  trueVal <- gamCoef[gamLoc]
  fx <- dgammaPost(gamVal=seq(from=0.05, to=5.00, by=0.05), gamCoef=gamCoef, 
                   gamLoc=gamLoc, W=W, X=X, U=U, xiDay=xiDay, hypGam=hypGam)
  width <- max(0.5, sum(fx > 0.01) * 0.05)
  domain <- c(max(0, trueVal - width), max(1, trueVal + width))
  
  plotGamMetrDen(gamMetrMat[,gamLoc], gamLoc, domain, numPts, trueVal, showDuns)
}




# Plots of the gamma posterior densities ---------------------------------------

numPts <- 100  # how many points to put in plots

plotGamPost(1)
plotGamPost(2)
plotGamPost(3)
plotGamPost(4)
plotGamPost(5)
plotGamPost(6)
plotGamPost(7)
plotGamPost(8)
plotGamPost(9)
plotGamPost(10)
plotGamPost(11)





