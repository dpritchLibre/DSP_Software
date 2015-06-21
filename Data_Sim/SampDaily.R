
sampDaily <- function(nTot, fwLen, varDist, varNames, fwDayNames=NULL) {
  p <- length(varDist)
  numDays <- nTot * fwLen
  if (is.null(fwDayNames))
    fwDayNames <- rep(paste0("day", 1:fwLen), times=nTot)
  
  dailyData <- lapply(X=1:p, FUN=function(x) replicate(n=nTot, expr=eval(varDist[x])) )
  names(dailyData) <- varNames

  return ( data.frame(cycleDay=fwDayNames, dailyData) )
}
