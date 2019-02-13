
#### pairwise permutation test for the difference of differentiation.1 among population from different regions (aggregates)
## between within
pwag.permutest=function(f,ncode,group, permutations){

  pairwisepermat=function (dat,group, permutations,  weight.type = 1, strata = NULL, parallel = getOption("mc.cores"))
  {
    EPS <- sqrt(.Machine$double.eps)
    classmean <- function(ind, dmat, indls) {
      sapply(indls, function(x)
        mean(c(dmat[ind == x, ind == x]),
             na.rm = TRUE))
    }
    mrpp.perms <- function(ind, dmat, indls, w) {
      weighted.mean(classmean(ind, dmat, indls), w = w, na.rm = TRUE)
    }
    if (inherits(dat, "dist"))
      dmat <- dat
    else if ((is.matrix(dat) || is.data.frame(dat)) &&
             isSymmetric(unname(as.matrix(dat)))) {
      dmat <- dat
      attr(dmat, "method") <- "user supplied square matrix"
    }
    else dmat <- as.dist(dat)
    if (any(dmat < -sqrt(.Machine$double.eps)))
      stop("dissimilarities must be non-negative")
    distance <- attr(dmat, "method")
    dmat <- as.matrix(dmat)
    diag(dmat) <- NA
    N <- nrow(dmat)
    group <- factor(group)
    indls <- levels(group)
    ncl <- sapply(indls, function(x) sum(group == x))
    w <- switch(weight.type, ncl, ncl - 1, ncl * (ncl - 1)/2)
    classdel <- classmean(group, dmat, indls)
    names(classdel) <- names(ncl) <- indls
    del <- weighted.mean(classdel, w = w, na.rm = TRUE)
    E.del <- mean(dmat, na.rm = TRUE)
    ## 'Classification strength' if weight.type == 1
    ## Do not calculate classification strength because there is no
    ## significance test for it. Keep the item in reserve for
    ## possible later re-inclusion.
    CS <- NA
    requireNamespace("vegan")
    permute::how
    permute::getBlocks
    permute::`setBlocks<-`
    permute::shuffleSet
    getPermuteMatrix=function (perm, N, strata = NULL)
    {
      if (length(perm) == 1) {
        perm <- how(nperm = perm)
      }
      if (!missing(strata) && !is.null(strata)) {
        if (inherits(perm, "how") && is.null(getBlocks(perm)))
          setBlocks(perm) <- strata
      }
      if (inherits(perm, "how"))
        perm <- shuffleSet(N, control = perm)
      else {
        if (!is.integer(perm) && !all(perm == round(perm)))
          stop("permutation matrix must be strictly integers: use round()")
      }
      if (is.null(attr(perm, "control")))
        attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"),
                                                nperm = nrow(perm)), class = "how")
      perm
    }
    permutations <- getPermuteMatrix(permutations, N, strata = strata)
    if (ncol(permutations) != N)
      stop(gettextf("'permutations' have %d columns, but data have %d rows",
                    ncol(permutations), N))

    control <- attr(permutations, "control")
    if(nrow(permutations)) {
      perms <- apply(permutations, 1, function(indx) group[indx])
      permutations <- ncol(perms)

      ## Parallel processing
      if (is.null(parallel))
        parallel <- 1
      hasClus <- inherits(parallel, "cluster")
      if (hasClus || parallel > 1) {
        if(.Platform$OS.type == "unix" && !hasClus) {
          m.ds <- unlist(mclapply(1:permutations, function(i, ...)
            mrpp.perms(perms[,i], dmat, indls, w),
            mc.cores = parallel))
        } else {
          if (!hasClus) {
            parallel <- makeCluster(parallel)
          }
          m.ds <- parCapply(parallel, perms, function(x)
            mrpp.perms(x, dmat, indls, w))
          if (!hasClus)
            stopCluster(parallel)
        }
      } else {
        m.ds <- apply(perms, 2, function(x) mrpp.perms(x, dmat, indls, w))
      }
      p <- (1 + sum(del + EPS >= m.ds))/(permutations + 1)
      r2 <- 1 - del/E.del
    } else { # no permutations
      m.ds <- p <- r2 <- NA
      permutations <- 0
    }
    out <- list(call = match.call(), delta = del, E.delta = E.del, CS = CS,
                n = ncl, classdelta = classdel, Pvalue = p, A = r2,
                distance = distance, weight.type = weight.type,
                boot.deltas = m.ds, permutations = permutations,
                control = control)
    class(out) <- "mrpp"
    out
  }
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
  return(pairwiseDav=pairwiseDav)
}

        pdis=diff1dis(f,ncode)
        rownames(pdis)=group
        colnames(pdis)=group

       # subset(pdis, rownames(pdis) %in% as.character(g1[c(1,4)]))
        g1=unique(group)
        dat1=matrix(data = 0, nrow = length(g1), ncol = length(g1))
        rownames(dat1)=g1
        colnames(dat1)=g1
        group=as.factor(group)
        group0=list()
        group00=list()
        pdis0=list()
        rcfilter=list()
        dat0=list()
          for (i in seq_along(g1) ) {
            group0[[i]]=list()
            group00[[i]]=list()
            pdis0[[i]]=list()
            rcfilter[[i]]=list()
            dat0[[i]]=list()
            for(j in seq_along(g1)){
              if (i==j) {
         #      print(c(i,j))
                dat1[i,j]=0}
              else {
                group0[[i]][[j]]=list()
                pdis0[[i]][[j]]=list()
                rcfilter[[i]][[j]]=list()
                dat0[[i]][[j]]=list()
                group0[[i]][[j]]=as.factor(c(as.character(group[group==g1[i]]),as.character(group[group==g1[j]])))
                group00[[i]][[j]]=sort(group0[[i]][[j]], decreasing =FALSE)
                # print(c(i,j))   as.character(group[group==g1[j]])))
                rcfilter[[i]][[j]]=which(rownames(pdis) %in% as.character(g1[c(i,j)]))
                pdis0[[i]][[j]]=pdis[rcfilter[[i]][[j]],rcfilter[[i]][[j]]]
                dat0[[i]][[j]]=pairwisepermat(pdis0[[i]][[j]],group =group00[[i]][[j]],permutations)
               dat1[i, j]=pairwisepermat(pdis0[[i]][[j]],group =group00[[i]][[j]],permutations)$Pvalue
              }
            }
          }




  ## will issue error if only a single group, this is not for pairwise
 # mod.aov <- anova(x)
 # nobs <- length(x$distances) ## number of observations
 # mod <- lm(x$distances ~ x$group)
 # mod.Q <- mod$qr
 # p <- mod.Q$rank
 # resids <- qr.resid(mod.Q, x$distances)

  ## extract groups
 # group <- x$group
 # ## permutations is either a single number, a how() structure or a
  ## permutation matrix
 # permutations <- getPermuteMatrix(permutations, nobs)
 # nperm <- nrow(permutations)

  ## pairwise comparisons
#  if(pairwise) {
  #  combin <- combn(levels(group), 2) # unique pairings
 #   n.pairs <- ncol(combin)
 # }

 return(list(pwpermutest.detail=dat0,pwtest.aggregate=dat1))

}

