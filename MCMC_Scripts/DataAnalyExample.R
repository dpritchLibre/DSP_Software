
# Read in data -----------------------------------------------------------------
{
  if (Sys.info()[8] == "dpritch") {
    setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
    outPath <- "/home/dpritch/Documents/Projects/Dunson Day Specific/Software/Output"
  }
  else {
    setwd("/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/")
    outPath <- "/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/Output"
  }
}

source("Format_Data/FormatDspDat.R")
source("MCMC_Scripts/McmcSampler.R")
load("Data/RawData.RData")




# Prepare data for use in MCMC sampler -----------------------------------------

idName <- "subjId"
cycName <- "cycle"
sexName <- "intercourse"
fwName <- "fwInd"
formula <- formula(pregInd ~ 0 + fwInd + age + bmi + gravid + cycleLen + lube)
fwLen <- 5

# Delete some arbitrary data to make sure that program can put the data together properly
dspDat <- dspDat(formula, baseline[-c(1,2), ], cycle[-(20:25), ], daily[-tail(1:nrow(daily), n=30), ], 
                 idName, cycName, sexName, fwName, fwLen)
summary(dspDat)

dspOut <- dsp(dspDat, numSamp=5e3)







# Summary stats of the output --------------------------------------------------

# Gamma coefs
# gamTab <- apply(read.csv(file=paste0(outPath, "GAMMA.csv")), MARGIN=2, 
#                 FUN=quantile, probs=c(0.025, 0.500, 0.975))
gamTab <- apply(dspOut$gam, MARGIN=2, FUN=quantile, probs=c(0.025, 0.500, 0.975))
trueVals <- c(0.14, 0.08, 0.34, 0.31, 0.08, exp(c(-0.08, -0.43, -1.03, -0.22, -0.47, 1.21)), rep(1,3))
rbind(gamTab, trueVals)

