
sampDataComb <- function(baselinePre, cyclePre, dailyPre) {
#   subjDays <- table(cycle$subjId) * fwLen  
# 
#   baseline <- data.frame( sapply(X=baseline, FUN=rep, times=subjDays, simplify=FALSE) )
#   cycle <- apply(X=cycle, MARGIN=2, FUN=rep, each=fwLen)
# 
#   return ( data.frame(cycle, baseline, daily) )
  cycLen <- Filter(dailyPre$cycLen, f=function(x) !is.na(x))
  basToCycExpan <- rep(1:nrow(baselinePre), times=table(cyclePre$subjId))
  basExpan <- rep(basToCycExpan, times=cycLen)
  cycExpan <- rep(1:nrow(cyclePre), times=cycLen)
  
  return ( data.frame( baselinePre[basExpan, ],
                       cyclePre[cycExpan, ],
                       dailyPre ) )
}


