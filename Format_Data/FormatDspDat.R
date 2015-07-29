
# Mung data into format for use by mcmc sampler ================================

dspDat <- function(baseline=NULL, cycle=NULL, daily, idName, cycName, 
                   pregName, sexName, fwIndName=NULL, varInclNames, fwLen) {
  source("Format_Data/FormatCheckValidInput.R")
  source("Format_Data/FormatHelperFcns.R")
  source("Format_Data/FormatCleanDat.R")
  source("Format_Data/FormatGetRedDat.R")
  source("Format_Data/FormatGetModelObj.R")
  source("Format_Data/FormatGetSamplerObj.R")
  source("Format_Data/FormatGetDatInfo.R")
  
  # Create vectors of var names needed for each dataset
  keepIfPregInSet <- function(set) if (pregName %in% names(set)) pregName else NULL
  varNames <- list( id = idName,
                    cyc = cycName,
                    preg = pregName,
                    sex = sexName,
                    fwInd = fwIndName,
                    basIncl = c(idName, varInclNames$baseline),
                    cycIncl = c(idName, cycName, keepIfPregInSet(cycle), varInclNames$cycle),
                    dayIncl = c(idName, cycName, keepIfPregInSet(daily), 
                                sexName, fwIndName, varInclNames$daily) )
  
  # TODO: sort data by id/cyc
  # TODO: check valid input
  
  # For daily data: remove non-FW days, cycles that have wrong number of FW days or include
  # missing in the cycle (in the model variables).  For baseline / cycle: remove observations 
  # that have missing data (in the model variables).
  cleanDat <- getCleanDat(baseline, cycle, daily, varNames, fwLen, cycInDaily)
  
  # Reduce data to subjects and cycles that are common to all datasets
  idVec <- getCommonId(cleanDat, idName)
  cycList <- getCommonCyc(cleanDat, varNames, idVec)
  redDat <- getRedDat(cleanDat, varNames, idVec, cycList)

  # Create X, Y, and U (from the Dunson and Stanford paper)
  modelObj <- getModelObj(redDat, varNames, fwLen, cycList)
  
  # Create objects for use in MCMC sampler (see 'FormatGetSamplerObj.R' for more details)
  samplerObj <- getSamplerObj(modelObj, fwLen)
  
  # Stats related to munging process for use by summary fcn
  datInfo <- getDatInfo(baseline, cycle, daily, cleanDat, redDat, 
                        modelObj, varNames, varInclNames, fwLen, idVec, cycList)
  
  # Construct dspDat object
  dspDat <- list( baseline = baseline,
                  cycle = cycle,
                  daily = daily,
                  cleanDat = cleanDat,
                  redDat = redDat,
                  modelObj = modelObj,
                  samplerObj = samplerObj,
                  datInfo = datInfo )
  structure(dspDat, class="dspDat")
}




# Describes the data munging process ===========================================

summary.dspDat <- function(dspDat) {
  datInfo <- dspDat$datInfo
  hline <- paste0(rep("-", 60), collapse="")
  
  # internal functions to assist printing --------------------------------------
  
  printStats <- function(dat) {
    cat("\n", hline, "\nRaw data:\n\n",
        "    baseline data:  ", dat$bas$sub, " subjects\n",
        "    cycle data:     ", dat$cyc$sub, " subjects with ", dat$cyc$cyc, " cycles\n",
        "    daily data:     ", dat$day$sub, " subjects with ", dat$day$cyc, " cycles and ", 
        dat$day$day, " days\n", sep="")
  }
  
  printVars <- function(charVec, returnWidth=45) {
    if (is.null(charVec))
      return (NULL)
    
    currLen <- 0
    outVec <- NULL
    
    for (i in 1:length(charVec)) {
      if (currLen >= returnWidth) {
        outVec <- paste0(outVec, "\n", paste(rep("", 22), collapse=" "))
        currLen <- 0
      }
      else if (currLen != 0) {
        outVec <- paste0(outVec, ", ", charVec[i])
        currLen <- currLen + nchar(charVec[i])
      }
      else {
        outVec <- paste0(outVec, charVec[i])
        currLen <- currLen + nchar(charVec[i])
      }
    }
    
    return (outVec)
  }
  
  # ----------------------------------------------------------------------------
  
  printStats(datInfo$numRaw)
  printStats(datInfo$numClean)

  numRed <- datInfo$numRed
  cat("\n", hline, "\nCombining the data:\n\n",
      "    common observations:  ", numRed$sub, 
      " subjects with ", numRed$cyc, 
      " cycles and ", numRed$day, " days\n", sep="")
  
  cat("\n", hline, "\nThe model variables:\n\n",
      "    variable names:  ", printVars(datInfo$modelVars), "\n",
      "    design matrix:   ", printVars(datInfo$designMatVars), "\n", sep="")
  
  cat(hline, "\n", sep="")
}


