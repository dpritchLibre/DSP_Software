
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
  # checkValidInput()
  
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
  
  numClean <- datInfo$numClean
  cat(hline, "\nAfter cleaning the data:\n\n",
      "    baseline data:  ", numClean$id$bas, " subjects\n",
      "    cycle data:     ", numClean$id$cyc, " subjects with ", numClean$cyc$cyc, " cycles\n",
      "    daily data:     ", numClean$id$day, " subjects with ", numClean$cyc$day, " cycles\n\n",
      sep="")
  
  numRed <- datInfo$numRed
  cat(hline, "\nCombining the data:\n\n",
      "    common observations:  ", numRed$subj, 
      " subjects with ", numRed$cyc, " cycles\n\n", sep="")
}


