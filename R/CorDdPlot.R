
CorDdPlot = function(x, d, ncode) {
    data = COR_DeltaDd(x, d, ncode)
   # requireNamespace("ggplot2")
     a= plot(data$Dgeo, data$PairwiseDeltaD, xlab = "Geographic Distance", ylab = "Genetic differentiation (DeltaD)")
    abline(lm(data$PairwiseDeltaD ~ data$Dgeo))
    lm = lm(data$PairwiseDeltaD ~ data$Dgeo)
    return(summary(lm))
}
