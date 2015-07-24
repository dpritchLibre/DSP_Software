
# Mung data into format for use by mcmc sampler ================================

dspDat <- function(baseline=NULL, cycle=NULL, daily, idName, cycleName, cycleDayName, 
                   pregName, sexName, varInclNames, fwDays) {
  fwLen <- length(fwDays)
  # TODO: check valid input
  cleanDat <- list(baselineCl=baseline[-(1:2),], cycleCl=cycle, dailyCl=daily)
  
  idVec <- getCommonId(cleanDat, idName)
  cycList <- getCommonCyc(cleanDat, idName, cycleName, idVec)
  redDat <- getRedDat(cleanDat, idName, cycleName, cycleDayName, pregName, sexName, 
                      varInclNames, idVec, cycList)

  modelObj <- getModelObj(redDat, idName, cycleName, cycleDayName, 
                          pregName, sexName, varInclNames, fwLen, cycList)
  samplerObj <- getSamplerObj(modelObj, fwLen)
  
  datInfo <- getDatInfo(cleanDat, redDat, idName, idVec, cycList, fwLen)
  dspDat <- list( cleanDat = cleanDat,
                  redDat = redDat,
                  modelObj = modelObj,
                  samplerObj = samplerObj,
                  datInfo = datInfo )
  
  structure(dspDat, class="dspDat")
}




# Describes the data munging process -------------------------------------------

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




# Check if input is valid ------------------------------------------------------
#
# no copies of id in baseline
# no copies of id+cyc in cycle
# no copies of id+cyc+day in daily
#
# baseline/cycle are NULL/mat/arr/df
# daily is mat/arr/df'
#
# X, Y 0/1 or yes/no
#
# (baseline == NULL) <==> (varInclNames$baseline == NULL and similarly for cycle
# varInclNames all match a name in corresp dataset and similarly for other *Name
# no multiple pregnancies
# no multiple id/cyc in cycle or id/cyc/day in day
# if preg in daily then preg consistent throughout cycle
#
# cycle and cycleDay are numeric
# no duplicates in fwDays




# Obtain vector of common id's that are needed ---------------------------------
#
# NULL[[idName]] == NULL
# unique(NULL) == NULL
# length(NULL) == 0

getCommonId <- function(cleanDat, idName) {
  list2env(cleanDat, envir=environment())
  basId <- baselineCl[[idName]]
  cycId <- unique( cycleCl[[idName]] )
  dayId <- unique( dailyCl[[idName]] )
  
  intrAllowNull <- function(a, b) if (!is.null(b)) intersect(a, b) else a  
  idVec <- Reduce(f=intrAllowNull, x=list(dayId, cycId, basId))
  
  return (idVec)
}




# Obtain list of common cycles that are needed ---------------------------------

getCommonCyc <- function(cleanDat, idName, cycleName, idVec) {
  list2env(cleanDat, envir=environment())
  n <- length(idVec)
  
  if (is.null(cycleCl))
    cycList <- lapply(idVec, function(x) unique(dailyCl[x == dailyCl[[idName]], cycleName]))
  else {
    cycById <- list( cyc = lapply(idVec, function(x) cycleCl[x == cycleCl[[idName]], cycleName]),
                     day = lapply(idVec, function(x) dailyCl[x == dailyCl[[idName]], cycleName]) )
    cycList <- lapply(1:n, function(j) intersect(cycById$cyc[[j]], cycById$day[[j]]) )
  }
  return (cycList)
}




# Reduce data sets to pertinant observations -----------------------------------

getRedDat <- function(cleanDat, idName, cycleName, cycleDayName, 
                      pregName, sexName, varInclNames, idVec, cycList) {
  list2env(cleanDat, envir=environment())
  pregInCycBool <- pregName %in% names(cycleCl)
  
  # Create reduced baseline data
  if (!is.null(baselineCl)) {
    baseBool <- (baselineCl[[idName]] %in% idVec)
    baselineRed <- baselineCl[baseBool, varInclNames$baseline, drop=FALSE]
  }
  else
    baselineRed <- NULL
  
  # Create reduced cycle data
  if (!is.null(cycleCl)) {
    cycObsIdx <- getObsIdx(cycleCl[[idName]], cycleCl[[cycleName]], idVec, cycList)
    cycKeepNames <- c(varInclNames$cycle, if (pregInCycBool) pregName else NULL)
    cycleRed <- cycleCl[cycObsIdx, cycKeepNames, drop=FALSE]
  }
  else
    cycleRed <- NULL

  # Create reduced daily data
  dayObsIdx <- getObsIdx(dailyCl[[idName]], dailyCl[[cycleName]], idVec, cycList)
  dayKeepNames <- c(idName, cycleName, cycleDayName, if (!pregInCycBool) pregName else NULL, 
                    sexName, varInclNames$daily)
  dailyRed <- dailyCl[dayObsIdx, dayKeepNames]
  
  return ( list( baselineRed = baselineRed, 
                 cycleRed = cycleRed, 
                 dailyRed = dailyRed ) )
}




# Obtain idx of common id's + cycles -------------------------------------------

getObsIdx <- function(datId, datCyc, idVec, cycList) {
  n <- length(idVec)
  idIdx <- lapply(idVec, function(x) which(datId == x))
  
  getCycBool <- function(j) ( datCyc[idIdx[[j]]] %in% cycList[[j]] )
  obsIdx <- lapply(1:n, function(j) subset(idIdx[[j]], getCycBool(j)) )
  
  return ( unlist(obsIdx) )
}




# Obtain info about the data munging process -----------------------------------

getDatInfo <- function(cleanDat, redDat, idName, idVec, cycList, fwLen) {
  list2env(cleanDat, envir=environment())
  list2env(redDat, envir=environment())
  
  numClean <- list( id = list( bas = length(baselineCl[[idName]]),
                               cyc = length( unique( cycleCl[[idName]] ) ),
                               day = length( unique( dailyCl[[idName]] ) ),
                               common = length(idVec) ),
                    cyc = list( cyc = nrow(cycleCl),
                                day = nrow(dailyCl) / fwLen ) )
  
  numRed <- list( subj = length(idVec),
                  cyc  = nrow(dailyRed) / fwLen )
  
  return ( list(numClean = numClean,
                numRed = numRed) )
}




# Combine base, cyc, day data into a daily dataset -----------------------------

getModelObj <- function(redDat, idName, cycleName, cycleDayName, 
                        pregName, sexName, varInclNames, fwLen, cycList) {
  list2env(redDat, envir=environment())
  n <- length( unique( dailyRed[[idName]] ) )
  N <- length( dailyRed[[idName]] )
  baseExpan <- rep(1:n, times=(sapply(X=cycList, FUN=length) * fwLen))

  # Combine datasets into a daily dataset that still contains factors
  covDat <- list( baselineRed[baseExpan, , drop=FALSE],
                  cycleRed[rep(1:(N / 5), each=fwLen), , drop=FALSE],
                  dailyRed[, c(cycleDayName, varInclNames$daily), drop=FALSE] )
  uFactor <- suppressWarnings( data.frame( Filter(length, covDat) ) ) # Row names conflict
  
  # Formula for use by model.matrix to convert factors to design matrix
  theModelFormula <- formula( paste0("~ -1 + ", cycleDayName, " + ", 
                                     paste(unlist(varInclNames), collapse=" + ")) )
  # Covariate matrix converted to design matrix
  U <- model.matrix(theModelFormula, data=uFactor)
  
  if (pregName %in% names(cycleRed))
    Y <- cycleRed[[pregName]]
  else
    Y <- dailyRed[seq(from=fwLen, to=N, by=fwLen), pregName]

  modelObj <- list( Y = Y,
                    X = dailyRed[[sexName]],
                    id = dailyRed[[idName]],
                    U = U )
  return (modelObj)
}




# Create objects for use in the MCMC sampler -----------------------------------

getSamplerObj <- function(modelObj, fwLen) {
  # Contains 'Y', 'X', 'id', 'U'
  list2env(modelObj, envir=environment())
  
  pregDayBool <- rep(convToBool(Y), each=fwLen)
  sexBool <- (convToBool(X) == 1)
  sexPregBool <- (sexBool & pregDayBool)
  
  n <- length(unique(id[sexBool])) # number of individuals
  q <- ncol(U) # number of covariates
  
  # Reduce data to days in which intercourse occured -----------------------------
  
  # Correspond to the rows after reducing data via sexBool / sexPregBool
  subSexRows <- replace(rep(0, length(id)), list=which(sexBool), values=1:sum(sexBool))
  subSexPregRows <- replace(rep(0, length(id)), list=which(sexPregBool), values=1:sum(sexPregBool))
  
  # Elements are indices corresponding to a subject
  idIdx <- Filter(length, lapply(X=unique(id), FUN=function(x) subSexRows[(id == x) & sexBool]))
  idDayExpan <- rep(1:length(idIdx), times=sapply(idIdx, length))
  
  # Elements are indices of cycles that have pregnancy
  pregCycIdx <- Filter(length, lapply(X=unique(id), FUN=function(x) 
    subSexPregRows[(id == x) & sexBool & pregDayBool]))

  # Reduce objects to intercourse days
  U <- U[sexBool, ]
  pregDayBool <- pregDayBool[sexBool]
  idPregIdx <- which( tapply(pregDayBool, INDEX=id[sexBool], FUN=function(x) TRUE %in% x) )

  # Convert binary cols of U to boolean
  covIsBinVec <- apply(U, MARGIN=2, FUN=function(x) all.equal(names(table(x)), c("0","1")))
  uBool <- lapply(1:q, function(j) if (!covIsBinVec[j]) NULL else (U[, j] == 1))
  pregUBool <- lapply(uBool, function(x) x[pregDayBool])
  
  # ----------------------------------------------------------------------------
    
  samplerObj <- list( U = U,
                      pregDayBool = pregDayBool,
                      pregCycIdx = pregCycIdx,
                      uBool = uBool,
                      pregUBool = pregUBool,
                      idIdx = idIdx,
                      idDayExpan = idDayExpan,
                      idPregIdx = idPregIdx,
                      n = n,
                      q = q )
  return (samplerObj)
}




# Convert possible 0/1 or "yes" / "no" to boolean ------------------------------
#
# PRE: already checked that all values are logical or 0/1 or start with "y","Y","n", or "N"

convToBool <- function(x) {
  if (is.logical(x))
    return (x)
  else if (is.numeric(x))
    return (x == 1)
  else
    sapply(substr(x, 1, 1), function(let) if ((let == "n") | (let == "N")) FALSE else TRUE)
}



