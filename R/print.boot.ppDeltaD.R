
print_boot_ppDeltaD<-function(x,...){
  a<-dim(x$lowerl)[1]
  ci<-matrix(nrow=a,ncol=a)
  dimnames(ci)<-dimnames(x$lowerl)
  for (i in 2:a)
    for (j in 1:(a-1)){
      ci[j,i]<-x$uperl[j,i]
      ci[i,j]<-x$lowerl[j,i]
    }
  cat("\n       Upper limit above diagonal \n")
  cat("       Lower limit below diagonal \n \n")
  print(ci)
  invisible(x)
}
