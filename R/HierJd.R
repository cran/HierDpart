#### This is the script for calculating the dissimilarity based on 1H, richness, which is the proportation of non-shared alleles

### This is used for the beta measures of 1H, compared to Fst, Delta D

HierJd=function(f,ncode,nreg, r){

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
    adegenet::locNames
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

  genfiles = read.genepop1(f, ncode=ncode, quiet = TRUE)  # covert the genepop #files to genind files, we can also use read.genpop from adegent package

  pairwise.propShared <- function(gi)
  {
    n.pops <- length(unique(pop(gi)))
    allPairs <- combn(1:n.pops, 2)
    gen.mat<- matrix(0, nrow=n.pops, ncol=n.pops)
    pops <- seppop(gi)
    pspop <- function(x)
    {
      pp  <- seppop(x)
      p1 <- pp[[1]]
      p2 <- pp[[2]]

      na <- ncol(p1@tab)
      maf <- NA
      m1 <- colMeans(p1@tab[,], na.rm=T)/2
      m2 <- colMeans(p2@tab[,], na.rm=T)/2

      m12 <- apply(rbind(m1,m2), 2, min, na.rm=T)

      lfl <- NA
      facs <- levels(p1@loc.fac)
      for (i in 1:length(locNames(p1))) 	lfl[i] <- sum(m12[p1@loc.fac==facs[i]])
      mean(lfl, na.rm=T)
    }

    for (i in 1:ncol(allPairs))
    {
      np1 <- allPairs[1,i]
      np2 <- allPairs[2,i]

      p12 <- repool(pops[[np1]], pops[[np2]])
      ps <- pspop(p12)
      gen.mat[np1,np2] <- ps
      gen.mat[np2,np1] <- ps

    }
    la <- levels(pop(gi))
    colnames(gen.mat) <- rownames(gen.mat) <- la
    return(as.dist(gen.mat))
  }

  Jdpop=pairwise.propShared(genfiles)
  Jd_pop_region=mean(Jdpop)

 # requireNamespace("dplyr")
  npops = length(levels(genfiles$pop))
  nloci = length(levels(genfiles$loc.fac))
  sampsize = summary(genfiles$pop)  ## sample size
  if (length(r) != nreg)
    stop("Number of regions should be equal to the number defined in the level")  ## number of pops per region
  if (sum(r) != npops)
    stop("Number of pops should be equal to the number defined in level")
  rsample = list()
  for (i in 1:nreg) {
    rsample[[i]] = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r, i)))])
  }
  rsample = as.data.frame(rsample)
  rsample = as.numeric(unlist(rsample))

  popr=list()
  for (i in 1:nreg) {
    popr[[i]] = list()
    popr[[i]] = as.factor(rep(paste("popr", i), times = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r, i)))])))  ### be aware that times depend on the sample size and str on your data
   # rsample[[i]] = sum(sampsize[(sum(head(r, i - 1)) + 1):(sum(head(r, i)))])
  }
  pop_region = unlist(popr)
  rgenfiles=genfiles
  rgenfiles$pop=pop_region
  Jdr=pairwise.propShared(rgenfiles)
  Jd_region_total=mean(Jdr)
  HierJd=cbind(Jd_region_total,Jd_pop_region)

  return(list(Jdpop=Jdpop,Jdr=Jdr,HierJd=HierJd))
}
