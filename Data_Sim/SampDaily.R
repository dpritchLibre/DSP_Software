
sampDaily <- function(nTot, varDist, varNames, cycLenDist, fwRevDays) {
  p <- length(varDist)
  
  # Sample the length of each cycle
  cycLen <- replicate(n=nTot, expr=eval(cycLenDist))
  
  # Sample the daily data covariates
  dailyData <- lapply(1:p, function(j) replicate(n=sum(cycLen), expr=eval(varDist[j])) )
  names(dailyData) <- varNames

  # Obtain the FW days by (i) calculating a reverse day number (by counting backwards from the 
  # last day of the cycle as day -1, -2, -3, etc), and (ii) then matching the reverse day to the
  # 'fwRevDays' chosen by the user
  cycleDay <- unlist( lapply(cycLen, FUN=seq, from=1) )
  reverseDay <- cycleDay - rep(cycLen + 1, times=cycLen)
  fwInd <- (reverseDay %in% fwRevDays) + 0
  
  # 'cycNum' is a vector of all NA's except for the last day of each cycle which contains
  # the value of the number of days in the cycle.  Thus the number of days in each cycle
  # can easily be extracted by removing the NA's from the vector
  cycLenIdx <- sapply(1:nTot, function(j) sum(cycLen[1:j]))
  cycLenInDaily <- replace(rep(NA, nTot), list=cycLenIdx, values=cycLen)
  
  return ( data.frame(dailyData, cycLen=cycLenInDaily, 
                      cycleDay=cycleDay, reverseDay=reverseDay, fwInd=fwInd) )
}
