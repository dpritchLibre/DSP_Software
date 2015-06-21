



sampDataComb <- function(baseline, cycle, daily, fwLen) {
  subjDays <- table(cycle$subjId) * fwLen  

  baseline <- data.frame( sapply(X=baseline, FUN=rep, times=subjDays, simplify=FALSE) )
  cycle <- apply(X=cycle, MARGIN=2, FUN=rep, each=fwLen)

  return ( data.frame(cycle, baseline, daily) )
}


