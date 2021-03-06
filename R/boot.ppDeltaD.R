## bootstraping of pairwise allelic differentiation
boot.ppDeltaD=function(dat,ncode,nboot,quant=c(0.025,0.975),diploid=TRUE,...){
  cl<-match.call()
 read.genepop1 <- function(file, ncode, quiet = TRUE) {
      adegenet::.readExt
      adegenet::.genlab
      adegenet::df2genind
      adegenet::is.genind
      adegenet::pop
      adegenet::repool
      adegenet::Hs
      adegenet::seppop
      adegenet::popNames
        if (toupper(.readExt(file)) != "GEN")
            stop("File extension .gen expected")
        if (!quiet)
            cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")
        prevcall <- match.call()
        txt <- scan(file, sep = "\n", what = "character", quiet = TRUE)
        if (!quiet)
            cat("\nFile description: ", txt[1], "\n")
        txt <- txt[-1]
        txt <- gsub("\t", " ", txt)
        NA.char <- paste(rep("0", ncode), collapse = "")
        locinfo.idx <- 1:(min(grep("POP", toupper(txt))) - 1)
        locinfo <- txt[locinfo.idx]
        locinfo <- paste(locinfo, collapse = ",")
        loc.names <- unlist(strsplit(locinfo, "([,]|[\n])+"))
        loc.names <- trimws(loc.names)
        nloc <- length(loc.names)
        txt <- txt[-locinfo.idx]
        pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
        npop <- length(pop.idx)
        nocomma <- which(!(1:length(txt)) %in% grep(",", txt))
        splited <- nocomma[which(!nocomma %in% pop.idx)]
        if (length(splited) > 0) {
            for (i in sort(splited, decreasing = TRUE)) {
                txt[i - 1] <- paste(txt[i - 1], txt[i], sep = " ")
            }
            txt <- txt[-splited]
        }
        pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))
        txt[length(txt) + 1] <- "POP"
        nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$", toupper(txt))) - 1
        pop <- factor(rep(1:npop, nind.bypop))
        txt <- txt[-c(pop.idx, length(txt))]
        temp <- sapply(1:length(txt), function(i) strsplit(txt[i], ","))
        ind.names <- vapply(temp, function(e) e[1], character(1))
        ind.names <- trimws(ind.names)
        vec.genot <- vapply(temp, function(e) e[2], character(1))
        vec.genot <- trimws(vec.genot)
        X <- matrix(unlist(strsplit(vec.genot, "[[:space:]]+")), ncol = nloc, byrow = TRUE)
        if (any(duplicated(ind.names))) {
            rownames(X) <- .genlab("", nrow(X))
        } else {
            rownames(X) <- ind.names
        }
        colnames(X) <- loc.names
        pop.names.idx <- cumsum(table(pop))
        pop.names <- ind.names[pop.names.idx]
        levels(pop) <- pop.names
        if (!all(unique(nchar(X)) == (ncode * 2)))
            stop(paste("some alleles are not encoded with", ncode, "characters\nCheck 'ncode' argument"))
        res <- df2genind(X = X, pop = as.character(pop), ploidy = 2, ncode = ncode, NA.char = NA.char)
        res@call <- prevcall
        if (!quiet)
            cat("\n...done.\n\n")
        return(res)
    }

  dat=read.genepop1(dat, ncode)
  if (adegenet::is.genind(dat)) dat1<-genind2hierfstat(dat)
  typ<-dim(dat)
  if(length(dim(dat1))==2){
    #Pop<-dat[,1]
    npop<-length(levels(dat$pop))
    nloc<-length(levels(dat$loc.fac))
    #ppsl<-array(numeric(npop*npop*nloc*3),dim=c(npop,npop,nloc,3))
    x0<-unique(dat1[,1])
    if(is.factor(dat1[,1])) dat1[,1]<-as.integer(dat1[,1])


 allele.freq<-function(population){
      # package require adegenet and pegas
      if (class(population) != "genind") {
        stop("You did not provide a valid genind object! Script stopped!")
      }

      # initial steps...
      numloci<-length(levels(population$loc.fac))  # this gets the total number of loci across all pops
      numpops<-length(levels(population$pop)) # this gets the total number of pops
      popnumallele<-population@loc.n.all     # this is a list of the population wide number of alleles at each pop
      lociname<-attributes(popnumallele)[[1]] # this is a list of the locinames (just L01, L02, L03,...)
      subdivpops<-seppop(population)

      # create list of matrices in which to place the numbers from summary
      alleletable<-vector("list",numloci)
      fralleletable<-vector("list",numloci)
      for(i in 1:numloci){
        alleletable[[i]]<-matrix(nrow=popnumallele[[i]],ncol=numpops)
        colnames(alleletable[[i]])<-popNames(population)
        rownames(alleletable[[i]])<-population@all.names[[i]][order(population@all.names[[i]])]
        fralleletable[[i]]<-matrix(nrow=popnumallele[[i]],ncol=numpops)
        colnames(fralleletable[[i]])<-popNames(population)
        rownames(fralleletable[[i]])<-population@all.names[[i]][order(population@all.names[[i]])]
      }

      # this is going to loop over all populations
      pegas::as.loci
      for (i in 1:numpops){
        x<-as.loci(subdivpops[[i]])
        s<-summary(x)
        # this loops over the loci
        for (j in 1:numloci){
          # this is the number of
          namevec<-(names(s[[j]]$allele))
          numnames<-length(namevec)
          # j<-2
          tablenames<-(rownames(alleletable[[j]]))
          for (k in 1:numnames){
            rownum<-which(tablenames==namevec[k])
            #  message("i = ",i," j = ",j," k = ",k," rownum = ",rownum)
            alleletable[[j]][rownum,i]<-s[[j]]$allele[k]
          }
        }
      }

      allpops<-as.loci(population)
      numbers<-summary(allpops)
      checkcnts<-matrix(nrow=numloci,ncol=2)
      for (i in 1:numloci){
        checkcnts[i,1]<-sum(numbers[[i]]$allele)
        checkcnts[i,2]<-sum(alleletable[[i]],na.rm=TRUE)
      }

      for (i in 1:numloci){
        for (j in 1:numpops){
          colsum<-sum(alleletable[[i]][,j],na.rm=TRUE)
          fralleletable[[i]][,j]<-round(alleletable[[i]][,j]/colsum, digits=3)
        }
      }

      freq=rapply(fralleletable, function(x) ifelse(is.na(x), 0, x), how="replace")

      alleletables<-list(count=alleletable, frequency=fralleletable,freq=freq)
      return(alleletables)
    }

    diff1dis=function(f){

      fr=allele.freq(f)
      af = fr$freq

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

      npop = length(levels(f$pop))
      nloc = length(levels(f$loc.fac))
      ppD<-array(numeric(npop*npop*nloc),dim=c(npop,npop,nloc))
      for( l in 1:nloc){
        for (i in 1:(npop-1)) {
          for (j in (i+1):npop) {
            #cat(i," ",j,"\n") #for debugging
            ppD[i,j,l]<-DeltaD((af[[l]][, c(i, j)]), str)[2]
          }
        }
      }
      #Dmat = list()
      #  for (l in 1:nloci) {
      #    Dmat[[l]] = matrix(data = 0, nrow = npops, ncol = npops)
      #    for (i in 1:npops) {
      #     for (j in 1:npops) {
      #       Dmat[[l]][i,j] = DeltaD((af[[l]][, c(i, j)]), str)[2]  ### select two pops from allelefrequency
      #
      #       }
      #   }
      #}
      return(ppD)
    }
    ppDeltaD=diff1dis(dat)
  }
  else
  {
    npop<-typ[1]
    nloc<-typ[3]
    ppDeltaD<-dat1
  }

  bppDeltaD<-array(numeric(npop*npop*nboot),dim=c(npop,npop,nboot))
  for (i in 1:nboot){
    sloc<-sample(nloc,replace=TRUE)
    ## sample loci and sum all of the DeltaD
    ovel<-apply(ppDeltaD[,,sloc],c(1,2),sum)
    bppDeltaD[,,i]<-ovel[,]/length(sloc)
  }
  #browser()
  stats::quantile
  lowerl<-apply(bppDeltaD,c(1,2),quantile,quant[1],na.rm=TRUE)
  uperl<-apply(bppDeltaD,c(1,2),quantile,quant[2],na.rm=TRUE)
  dimnames(lowerl)[[1]]<-dimnames(lowerl)[[2]]<-sort(x0)
  dimnames(uperl)<-dimnames(lowerl)
  res<-(list(call=cl,lowerl=lowerl,uperl=uperl,DeltaD.per.loc=ppDeltaD))
  class(res)<-"boot.ppDeltaD"
  res
}
