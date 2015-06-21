


sampCycle <- function(subjId, entryDist, lengthDist, maxStudyCyc=NULL, varDist, varNames) {
  
  subjCyc <- sampSubjCyc(subjId, entryDist, lengthDist, maxStudyCyc)
  cycVars <- sampCycVars(nrow(subjCyc), varDist, varNames)
  
  return ( data.frame(subjCyc, cycVars) )
}


sampSubjCyc <- function(subjId, entryDist, lengthDist, maxStudyCyc=NULL) {
  n <- length(subjId)
  
  # TODO: validate distribution expressions
  
  cycStartVec <- replicate(n=n, expr=eval(entryDist))
    # cycleLenVec interpreted as *additional cycles in study after the first* as a convenience
    # so that sampling a 0 makes sense for distribs such as geometric or Poisson
  cycLenVec <- replicate(n=n, expr=eval(lengthDist))
  
  if (is.null(maxStudyCyc))
    cycEndVec <- cycStartVec + cycLenVec
  else
    cycEndVec <- pmin(cycStartVec + cycLenVec, maxStudyCyc)
  
  subjId <- rep.int(x=subjId, times=(cycEndVec - cycStartVec + 1))
  cycle <- unlist( mapply(FUN=seq.int, from=cycStartVec, to=cycEndVec) )
  
  return ( data.frame(subjId=subjId, cycle=cycle) )
}


sampCycVars <- function(nTot, varDist, varNames=NULL) {
  p <- length(varDist)
  
  cycData <- lapply(X=1:p, FUN=function(x) replicate(n=nTot, expr=eval( varDist[[ x ]] )))
  names(cycData) <- varNames
  
  return ( data.frame(cycData) )
}
                       