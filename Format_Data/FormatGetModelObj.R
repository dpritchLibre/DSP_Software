
# Combine base, cyc, day data into desing matrix -------------------------------

getModelObj <- function(redDat, varNames, varInclNames, fwLen, cycList) {
  niVec <- Filter(sapply(cycList, length), f=as.logical)
  numOf <- list( subj = length(niVec),
                 cyc  = sum(niVec) )
             
  basExpan <- rep(1:numOf$subj, times=(niVec * fwLen))
  cycExpan <- rep(1:numOf$cyc, each=fwLen)
  
  # Combine datasets into a daily dataset that still contains factors
  fwDay <- as.factor( rep(1:fwLen, times=numOf$cyc) )
  dayKeepVars <- setdiff(varNames$dayIncl, c(varNames$id, varNames$cyc))
  covDat <- list( redDat$bas[basExpan, , drop=FALSE],
                  redDat$cyc[cycExpan, , drop=FALSE],
                  redDat$day[, , drop=FALSE] )
  # Suppressed warning: rows have same names (due to expansion)
  uFactor <- suppressWarnings( data.frame(fwDay, Filter(length, covDat)) ) 
  
  # Formula for use by model.matrix to convert factors to design matrix
  theModelFormula <- formula( paste(c("~ -1 + fwDay", unlist(varInclNames)), collapse=" + ") )
  # Covariate matrix converted to design matrix
  U <- model.matrix(theModelFormula, data=uFactor)
  
  if (varNames$preg %in% names(redDat$cyc))
    Y <- redDat$cyc[[varNames$preg]]
  else
    Y <- redDat$day[seq(from=fwLen, to=nrow(redDat$day), by=fwLen), varNames$preg]
  
  modelObj <- list( Y = Y,
                    X = redDat$day[[varNames$sex]],
                    id = redDat$day[[varNames$id]],
                    U = U )
  return (modelObj)
}