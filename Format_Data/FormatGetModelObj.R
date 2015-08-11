
# Combine base, cyc, day data into desing matrix -------------------------------

getModelObj <- function(redDat, varNames, fwLen, cycList) {
  niVec <- Filter(sapply(cycList, length), f=as.logical)
  numOf <- list( subj = length(niVec),
                 cyc  = sum(niVec) )
             
  basExpan <- rep(1:numOf$subj, times=(niVec * fwLen))
  cycExpan <- rep(1:numOf$cyc, each=fwLen)
  
  # Combine datasets into a daily dataset that still contains factors
  redDat$day[[varNames$fw]] <- as.factor( rep(1:fwLen, times=numOf$cyc) )
  covDat <- Filter( length, list( redDat$bas[basExpan, , drop=FALSE],
                                  redDat$cyc[cycExpan, , drop=FALSE],
                                  redDat$day[, , drop=FALSE] ) )
  # Covariate matrix converted to design matrix
  U <- model.matrix(formula, data=data.frame(covDat))
  
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