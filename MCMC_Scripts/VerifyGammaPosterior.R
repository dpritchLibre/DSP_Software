
# Read in data -----------------------------------------------------------------
{
  if (Sys.info()[8] == "dpritch")
    setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
  else
    setwd("/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/")
}

load("Data/PracticeDatAuxMod.RData")
source("MCMC_Scripts/McmcGetLogR.R")
set.seed(0)




# Dataset attributes -----------------------------------------------------------

fwLen <- 5
U <- as.matrix( U[, 1:11] )  # Only keep terms with nonzero beta coefs



# Set (true) parameter values --------------------------------------------------

betaDays <- log( c(0.14, 0.08, 0.34, 0.31, 0.08) )
betaCovs <- list( age = c(-0.08, -0.43, -1.03),
                  bmi = c(-0.22, -0.47),
                  gravid = 1.21 )

betaCoef <- c(betaDays, unlist(betaCovs))
names(betaCoef)[1:5] <- paste0("day", 1:5)
phi <- 1
xiCyc <- rep(xi, times=(table(id) / fwLen))
xiDay <- rep(x=xiCyc, each=fwLen)




# Sample p(phi | ...) via Metropolis algorithm =================================

B <- 5e4
delta <- 0.2
currBeta <- log(0.34)
beta3 <- numeric(length=B)
acceptBool <- logical(length=B)
pJ <- 0.1
changeLoc <- 3
hyperparGam <- list(boundL=0, boundU=Inf, ptMassProp=0.5, shape=1, rate=1)


for (i in 1:B) {
  
  if (rbinom(n=1, size=1, prob=pJ) == 1)
    proposeBeta <- 0
  else 
    proposeBeta <- runif(n=1, min=(currBeta - delta), max=(currBeta + delta))
  
  logR <- getLogR(Y=Y, U=U, xi=xiDay, betaCoef=betaCoef, currVal=currBeta, proposeVal=proposeBeta, 
                  changeLoc=changeLoc, hyperparGam=hyperparGam)
  if (log(runif(1)) < logR) {
    beta3[i] <- currBeta <- proposeBeta
    acceptBool[i] <- TRUE
  }
}




# Calculate quantile function for p(phi | ...) =================================

gamLoc <- 3
gamCoef <- exp( betaCoef[-gamLoc] )
a <- 1                               # gamma_h shape hyperparam
b <- 1                               # gamma_h rate hyperparam
bndL <- 0                            # gamma_h lower bound
bndU <- Inf                          # gamma_h upper bound

x <- seq(from=0.01, to=1, by=0.01)
Fx <- numeric(length(x))
for (k in 1:length(x))
  Fx[k] <- pgammaPost(x[k], gamCoef, gamLoc, W, X, U, p, a, b, bndL, bndU)$value
Fx <- unlist(Fx)




# Compare Metropolis and theoretical -------------------------------------------


Fmet <- sapply(X=seq(0.01,1.00, by=0.01), FUN=function(x) mean(exp(beta3) <= x))
round(data.frame(Fx=Fx, Fmet=Fmet), digits=2)



