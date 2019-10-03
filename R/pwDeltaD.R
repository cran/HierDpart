#### pairwise Delta D for a data.frame

pwDeltaD=function(x,para){

  #x=NormalizeData(x)
  # x= as.data.frame(apply(x, 2, function(x) (x - min(x))/(max(x)-min(x))))

  x=t(as.matrix(x))
  DeltaD = function(abun, struc) {
    n = sum(abun,na.rm = TRUE)
    N = ncol(abun)
    ga = rowSums(abun,na.rm = TRUE)
    gp = ga[ga > 0]/n
    G = sum(-gp * log(gp))
    H = nrow(struc)
    A = numeric(H - 1)
    W = numeric(H - 1)
    Diff = numeric(H - 1)
    wi = colSums(abun,na.rm = TRUE)/n
    W[H - 1] = -sum(wi[wi > 0] * log(wi[wi > 0]),na.rm = TRUE)
    pi = sapply(1:N, function(k) abun[, k]/sum(abun[, k], na.rm = TRUE))
    Ai = sapply(1:N, function(k) -sum(pi[, k][pi[, k] > 0] * log(pi[, k][pi[, k] > 0]),na.rm = TRUE))
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
          }
          else {
            ai[, j] = rowSums(abun[, II],na.rm = TRUE)
          }
        }
        pi = sapply(1:NN, function(k) ai[, k]/sum(ai[,k],na.rm = TRUE))
        wi = colSums(ai)/sum(ai)
        W[i - 1] = -sum(wi * log(wi))
        Ai = sapply(1:NN, function(k) -sum(pi[, k][pi[,k] > 0] * log(pi[, k][pi[, k] > 0])))
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
  nvar = ncol(x)
  napops=colnames(x)
  # pairwise matrix index
  pw <- combn(nvar, 2,simplify = FALSE) #,simplify = FALSE

  # set up parallel cluster
  #if(para){

  #  ncor <- parallel::detectCores()
  #}
  #cl <- parallel::makeCluster(ncor)

  # Dmat = matrix(data = 0, nrow = nvar, ncol = nvar,dimnames = list(napops, napops))
  if(para){
    ncor <- parallel::detectCores()
    print(paste("Using parallel",ncor-2,"cores"))
    Dp=parallel::mclapply(pw,function(i, abun) DeltaD(abun[,i],struc=str)[2], abun = x,  mc.cores = 2)
  }
  else {
    print(paste("Using normal computer resource"))
    Dp=lapply(pw,function(i, abun) DeltaD(abun[,i],struc=str)[2], abun = x)
  }

  pwD <- diag(nvar)
  diag(pwD) <- 0
  pwD[lower.tri(pwD)] <- unlist(Dp)
  # this ugly looking bit is due to the matrix filling up by columns instead of rows
  pwD[upper.tri(pwD)] <- t(pwD)[upper.tri(t(pwD))]
  colnames(pwD)  <- colnames(x)
  row.names(pwD)=colnames(x)



  #  Dmat=lapply(x,DeltaD()[2] i,j)
  #  b=parallel::parLapply(cl,pw,function(i,abun) HierDpart::IDIP(abun[,i],struc=str), abun=dat)
  #x[i]

  # parallel::stopCluster(cl)

  ###  use this loop only for small data

  #   for (i in 2:nvar) {
  #   for (j in 1:(i-1)) {
  #    Dmat[i, j] =Dmat[j, i]=DeltaD((x[, c(i, j)]),str)[2]
  #}
  #  }
  #  Dm=parallel::parLapply(x,DeltaD(x,str)[2], i=2:nvar,j=1:(i-1))
  #parallel::stopCluster(cl)

  # pairwiseDav=as.matrix(Dmat)
  #colnames(pairwiseDav) = colnames(x)
  #rownames(pairwiseDav) = colnames(x)
  ##DeltaDmat = as.dist(pairwiseDav, diag = FALSE, upper = FALSE)
  # or may try sapply(1:N, function(x) DeltaD((x[, c(i, j)]), str)[2],)

  return(list(PairwiseDeltaD = pwD))
}

