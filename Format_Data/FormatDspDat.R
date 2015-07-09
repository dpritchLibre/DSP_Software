
# Mung data into format for use by mcmc sampler ================================

makeDspDat <- function(baseline, cycle, daily, idName, cycleName, 
                       cycleDayName, pregName, intercourseName, varInclNames, fwLen) {
  # TODO: sort data by id before processing?
  
  idVec <- getCommonId(baseline, cycle, daily, idName)
  cycleList <- getCommonCyc(cycle, daily, idName, cycleName, idVec)
    # Datasets reduced to obs / vars for use in model
  reducDatList <- getReducData(baseline, cycle, daily, idName, cycleName, cycleDayName, 
                               varInclNames, idVec, cycleList)
  
  dspDat <- makeMcmcDat(reducDatList, idName, cycleName, cycleDayName, 
                        intercourseName, varInclNames, fwLen, cycleList)
  return (dspDat)
}




# Obtain vector of common id's that are needed =================================

getCommonId <- function(baseline, cycle, daily, idName) {
  
  # TODO: write this function for when id's not identical across datasets
  return ( baseline[[idName]] )
}




# Obtain list of common cycles that are needed =================================

# TODO: assumes that data comes from both cycle and daily sets.  Write for when this 
# is not the case.  daily set already accounted for in idVec

getCommonCyc <- function(cycle, daily, idName, cycleName, idVec) {
  n <- length(idVec)
  
    # Lists iterating over subj's with vec's of observation idx's
  idIdxCyc <- lapply(X=1:n, FUN=function(x) which(idVec[x] == cycle[[idName]]))
  idIdxDay <- lapply(X=1:n, FUN=function(x) which(idVec[x] == daily[[idName]]))
    # Lists iterating over subj's with vec's of cycle numbers
  cycListCyc <- lapply(X=idIdxCyc, FUN=function(x) cycle[x,cycleName])
  cycListDay <- lapply(X=idIdxDay, FUN=function(x) unique(daily[x,cycleName]))

    # List iterating over subj's with cycles that are common to both cycle and daily sets
  cycleList <- lapply(X=1:n, FUN=function(x) intersect(cycListCyc[[x]], cycListDay[[x]]))
  return (cycleList)
}




# Reduce data sets to pertinant observations ===================================

getReducData <- function(baseline, cycle, daily, idName, cycleName, cycleDayName, 
                           varInclNames, idVec, cycleList) {
  
  cycBoolVec <- getObsBool(cycleVec=cycle[[idName]], obsVec=cycle[[cycleName]], 
                           idVec=idVec, cycleList=cycleList)
  dayBoolVec <- getObsBool(cycleVec=daily[[idName]], obsVec=daily[[cycleName]], 
                           idVec=idVec, cycleList=cycleList)
  
  # TODO: consider what to do if no vars from a dataset to be used in model
  baselineRed <- subset(x=baseline[baseline[[idName]] %in% idVec, ], select=varInclNames$baseline)
  cycleRed <- subset(x=cycle[cycBoolVec, ], select=varInclNames$cycle)
  dailyRed <- daily[dayBoolVec, c(idName, cycleName, cycleDayName, varInclNames$daily)]
  
  return ( list( baseRed = baselineRed, 
                 cycRed = cycleRed, 
                 dayRed = dailyRed ) )
}




# Obtain bool with rows of common ids / cycles =================================

getObsBool <- function(cycleVec, obsVec, idVec, cycleList) {
  n <- length(idVec)
  
    # Cols in matrix iterate over subject and return whether obs needed in reduced data for subj
  cycBoolMat <- sapply(X=1:n, FUN=function(x) 
    (cycleVec %in% idVec[x]) & (obsVec %in% cycleList[[x]]) )
  
    # Takes union across cols and hence obtains whether each obs needed over all subjects
  cycBool <- apply(X=cycBoolMat, MARGIN=1, FUN=function(x) TRUE %in% x)
  return (cycBool)
}




# Combine base, cyc, day data into a daily dataset =============================

makeMcmcDat <- function(reducDatList, idName, cycleName, cycleDayName, 
                        intercourseName, varInclNames, fwLen, cycleList) {
  n <- length( unique( daily[[idName]] ) )  
    # TODO: rewrite using list2env
  baseRed <- reducDatList$baseRed
  cycRed <- reducDatList$cycRed
  dayRed <- reducDatList$dayRed
  
  n <- length( baseline[[idName]] )
    # Index by which to expand daily into cycle dataset
  baseExpandIdx <- unlist( mapply(FUN=rep, x=1:n, each=sapply(X=cycleList, FUN=length) ))
    # Expand baseline and cycle datasets into a daily dataset  
  baseExpand <- lapply(X=baseRed[baseExpandIdx, ], FUN=rep, each=fwLen)
  cycExpand <- lapply(X=cycRed, FUN=rep, each=fwLen)

    # Combine datasets into a daily dataset that still contains factors
  uFactor <- data.frame( baseExpand,
                         cycExpand,
                         daily[c(cycleDayName, varInclNames$daily)] )
  
    # Formula for use by model.matrix to convert factors to design matrix
  theModelFormula <- formula( paste0("~ -1 + ", cycleDayName, " + ", 
                                    paste(unlist(varInclNames), collapse=" + ")) )
    # Covariate matrix converted to design matrix
  uDesignMat <- model.matrix(theModelFormula, data=uFactor)

  
  dspDat <- list( Y = cycle[[pregName]],
                  X = (daily[[intercourseName]] == "yes") + 0,
                  id = daily[[idName]],
                  U = uDesignMat )
  return (dspDat)
}



