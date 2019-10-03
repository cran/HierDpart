\name{pwDeltaD}
\alias{pwDeltaD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pairwise Shannon differentiation (Delta D)
%%  ~~function to do ... ~~
}
\description{ Estimating pairwise Shannon differentiation (Delta D)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pwDeltaD(x, para)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input file, data matrix or data frame
%%     ~~Describe \code{x} here~~
}
  \item{para}{ Indicating whether using parallel process
%%     ~~Describe \code{para} here~~
}
}
\details{ Calculating pairwise Shannon differentiation, Delta D (Gaggiotti, O. E.,et al, 2018).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{PairwiseDeltaD}{Pairwise Delta D}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary applications, 11(7), 1176-1193.
%% ~put references to the literature/web site here ~
}
\author{ qinxinghu@gmail.com
%%  ~~who you are~~
}
\note{ Note that for windows system, the parallel doesn't work for multiple cores.In that case, users should set para=FALSE.
%%  ~~further notes~~
}

\examples{
### Users can test the example in their own computer,
##as this is pairwise differentiation, using parallel
# will save more time
#f <- system.file('extdata',package='HierDpart')
#infile <- file.path(f, "Island.gen")
#library(adegenet)
#Islandfile=read.genepop(infile,ncode=3)
#Islandfreq=tab(Islandfile,freq=TRUE)
#IslandpwD=pwDeltaD(Islandfreq,para=FALSE)
#print(IslandpwD)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

