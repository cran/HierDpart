##### Plot the ordination of different aggregates

plot_pcoa_aggregate=function(f,ncode,level,label = TRUE) {

  diff1dis=function(f,ncode){
  diveRsity::readGenepop
  gp = ncode
  fr = readGenepop(f, gp, bootstrap = FALSE)
  af = fr$allele_freq

  DeltaD = function(abun, struc) {
    ## Chao et al, 2017
    n = sum(abun)
    N = ncol(abun)
    ga = rowSums(abun)
    gp = ga[ga > 0]/n
    G = sum(-gp * log(gp))
    H = nrow(struc)
    A = numeric(H - 1)
    W = numeric(H - 1)
    Diff = numeric(H - 1)
    wi = colSums(abun)/n
    W[H - 1] = -sum(wi[wi > 0] * log(wi[wi > 0]))
    pi = sapply(1:N, function(k) abun[, k]/sum(abun[, k]))
    Ai = sapply(1:N, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0])))
    A[H - 1] = sum(wi * Ai)
    if (H > 2) {
      for (i in 2:(H - 1)) {
        I = unique(struc[i, ])
        NN = length(I)
        ai = matrix(0, ncol = NN, nrow = nrow(abun))
        c
        for (j in 1:NN) {
          II = which(struc[i, ] == I[j])
          if (length(II) == 1) {
            ai[, j] = abun[, II]
          } else {
            ai[, j] = rowSums(abun[, II])
          }
        }
        pi = sapply(1:NN, function(k) ai[, k]/sum(ai[, k]))
        wi = colSums(ai)/sum(ai)
        W[i - 1] = -sum(wi * log(wi))
        Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0])))
        A[i - 1] = sum(wi * Ai)
      }
    }
    Diff[1] = (G - A[1])/W[1]
    if (H > 2) {
      for (i in 2:(H - 1)) {
        Diff[i] = (A[i - 1] - A[i])/(W[i] - W[i - 1])
      }
    }
    Diff = Diff
    out = matrix(c(Diff), ncol = 1)
    return(out)
  }

  v1 = c("ecosystem", "region1", "pop1")
  v2 = c("ecosystem", "region1", "pop2")
  str = data.frame(v1, v2)
  str = as.matrix(str)
  npops = fr$npops
  nloci = fr$nloci
  Dmat = list()
  for (l in 1:nloci) {
    Dmat[[l]] = matrix(data = 0, nrow = npops, ncol = npops)
    for (i in 1:npops) {
      for (j in 1:npops) {
        Dmat[[l]][i, j] = DeltaD((af[[l]][, c(i, j)]), str)[2]  ### select two pops from allelefrequency
      }
    }
  }
  pairwiseDav = Reduce("+", Dmat)/length(Dmat)
  colnames(pairwiseDav) = fr$pop_names
  rownames(pairwiseDav) = fr$pop_names
  # library(popbio)
  DeltaDmat = as.dist(pairwiseDav)
return(list(pairwiseDav=pairwiseDav,DeltaDmat=DeltaDmat))
}

diff1=diff1dis(f,ncode)
ape::pcoa
dd2=pcoa(diff1$pairwiseDav, correction="none", rn=rownames(diff1$pairwiseDav))

biplot(dd2,  plot.axes = c(1,2), col=1: length(levels(as.factor(level))), dir.axis1=1, dir.axis2=1, cex=NULL, main="Aggregation plot",xlabs = "PCoA1", ylabs = "PCoA2")
#ordihull(dd2$vectors, level, col=1:4, lwd=3)
ordiellipse(dd2$vectors, level, col=1: length(levels(as.factor(level))), kind = "ehull", lwd=3)
#ordiellipse(dd2$vectors, level, col=1:4, draw="polygon")
ordispider(dd2$vectors, groups =level, display = "sites", col=1: length(levels(as.factor(level))), label=label)
points(dd2$vectors, pch=21, col="purple", bg=c("white"), cex=1.3)

}

