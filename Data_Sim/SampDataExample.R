
# Load data simulation functions -----------------------------------------------
#
# TODO: encapsulate the files into 1 driver

{
if (Sys.info()[8] == "dpritch")
  setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
else
  setwd("sam's working directory please change")
}

source("Data_Sim/SampBaseline.R")
source("Data_Sim/SampCycle.R")
source("Data_Sim/SampDaily.R")
source("Data_Sim/SampDataComb.R")
source("Data_Sim/SampPreg.R")




# Define study id's (and hence size)  ------------------------------------------

subjId <- 1:300




# Sample baseline data ---------------------------------------------------------

baseNames <- c("age","race","bmi","gravid","edu","depr","smoke","drinkAlc","hormContrac")

baseDist <- expression( sample( c("< 30","30-35","35-38","38+"), size=1,
                                prob=c(0.41, 0.39, 0.11, 0.09) ),              #  <-- age
                        sample( c("cauc","black","hisp","other"), 
                                size=1, prob=c(0.64, 0.12, 0.16, 0.08) ),      #  <-- race
                        sample( c("< 25","25-30","30+"), size=1, 
                                prob=c(0.63, 0.21, 0.16) ),                    #  <-- bmi
                        sample( c("no","yes"), size=1, prob=c(0.42,0.58) ),    #  <-- gravid
                        sample( c("no college","college","grad/prof"), 
                                size=1, prob=c(0.58, 0.30, 0.12) ),            #  <-- edu
                        sample( c("no","yes"), size=1, prob=c(0.90, 0.10) ),   #  <-- depr
                        sample( c("no","yes"), size=1, prob=c(0.82, 0.18) ),   #  <-- smoke
                        sample( c("no","yes"), size=1, prob=c(0.44, 0.56) ),   #  <-- drinkAlc
                        sample( c("no","yes"), size=1, prob=c(0.83, 0.17) ) )  #  <-- hormContrac

baselinePre <- sampBaseline(subjId=subjId, varDist=baseDist, varNames=baseNames)




# Sample cycle data ------------------------------------------------------------

  # Describe the left and right truncation of observed data.  'lengthDist' is interpreted by
  # sampCycle as *additional cycles in study after the first* (i.e. 0 is an allowable value).

entryDist <- expression( sample.int(n=3, size=1, prob=c(0.34,0.33,0.33)) )
lengthDist <- expression( rgeom(n=1, prob=0.25) )


  # Describe the data-generating process for cycle-specific variables

cycNames <- c("cycleLen","opk_use","bleed_intermen","bleed_luteal")

cycDist <- expression( sample( c("< 28 days","28-31 days","32+ days"), 
                               size=1, prob=c(0.11, 0.69, 0.20) ),             #  <-- cycleLen
                       sample( c("no","yes"), size=1, prob=c(0.70, 0.30) ),    #  <-- opk_use
                       sample( c("no","yes"), size=1, prob=c(0.57, 0.43) ),    #  <-- bleed_intermen
                       sample( c("no","yes"), size=1, prob=c(0.64, 0.36) ) )   #  <-- bleed_luteal


  # Sample cycle data from the aforementioned distributions

cyclePre <- sampCycle(subjId=subjId, entryDist=entryDist, lengthDist=lengthDist,
                      maxStudyCyc=12, varDist=cycDist, varNames=cycNames)




# Sample daily data ------------------------------------------------------------

dailyNames <- c("intercourse","cm_monit","lube")

dailyDist <- expression( sample( c("no","yes"), size=1, prob=c(0.64, 0.36) ),    #  <-- intercourse
                         sample( c("didn't check","type 1",
                                   "type 2","type 3","type 4"),
                                 size=1, prob=c(0.04, 0.05, 0.08, 0.04, 0.79) ), #  <-- cm_monit
                         sample( c("no","yes"), size=1, prob=c(0.88, 0.12) ) )   #  <-- lube

dailyPre <- sampDaily(nTot=nrow(cyclePre), fwLen=5, varDist=dailyDist, varNames=dailyNames)




# Combine data into a daily set ------------------------------------------------

analyData <- sampDataComb(baseline=baselinePre, cycle=cyclePre, daily=dailyPre, fwLen=5)




# Sample subject pregnancy -----------------------------------------------------

betaDays <- log( c(0.14, 0.08, 0.34, 0.31, 0.08) )
betaCovs <- list( age = c(-0.08, -0.43, -1.03),
                  bmi = c(-0.22, -0.47),
                  gravid = 1.21 )

pregVecPre <- sampPreg(dspDat=analyData, betaDays=betaDays, betaCovs=betaCovs, phi=1, fwLen=5)

finalData <- rmSuperfluous(baselinePre, cyclePre, dailyPre, pregVecPre)
invisible( list2env(finalData, envir=.GlobalEnv) )
rm(list=setdiff(ls(), c("baseline","cycle","daily")))




# Format data into format for use by mcmc sampler ------------------------------

source("Format_Data/FormatDspDat.R")

idName <- "subjId"
cycleName <- "cycle"
cycleDayName <- "cycleDay"
pregName <- "pregInd"
intercourseName <- "intercourse"
varInclNames <- list( baseline = c("age","bmi","gravid"),
                      cycle = "cycleLen",
                      daily = "lube" )

dspDat <- makeDspDat(baseline, cycle, daily, idName, cycleName, cycleDayName, 
                     pregName, intercourseName, varInclNames, fwLen=5)
invisible( list2env(dspDat, envir=.GlobalEnv) )
save(Y, X, id, U, file="Data/PracticeDat.RData")






