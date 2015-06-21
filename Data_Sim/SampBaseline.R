

sampBaseline <- function(subjId, varDist, varNames=NULL) {
  
  ## TODO: Need to test if valid input ##
  
  n <- length(subjId)
  p <- length(varDist)

  if (is.null(varNames))
    varNames <- sapply(X=1:p, FUN=function(x) paste0("base", x))

  baseData <- lapply(X=1:p, FUN=function(x) replicate(n=n, expr=eval( varDist[[ x ]] )))
  names(baseData) <- varNames

  return ( data.frame(baseData) )
}