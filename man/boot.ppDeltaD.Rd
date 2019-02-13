\name{boot.ppDeltaD}
\alias{boot.ppDeltaD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bootstrapping test of pairwise allelic differentiation (Delta D) between populations
%%  ~~function to do ... ~~
}
\description{ Do bootstrapping on pairwise allelic differentiation (Delta D) between populations and get the confidence interval.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
boot.ppDeltaD(dat, ncode, nboot, quant = c(0.025, 0.975), diploid = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{ Genetic file
%%     ~~Describe \code{dat} here~~
}
  \item{ncode}{ Genotype coding type
%%     ~~Describe \code{ncode} here~~
}
  \item{nboot}{ Number of bootstraps
%%     ~~Describe \code{nboot} here~~
}
  \item{quant}{ confidence interval
%%     ~~Describe \code{quant} here~~
}
  \item{diploid}{Ploidy
%%     ~~Describe \code{diploid} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{Do bootstrapping on pairwise allelic differentiation (Delta D) between populations and get the confidence interval
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{lowerl }{lower confidence interval of pairwise allelic differentiation}
\item{uperl }{uper confidence interval of pairwise allelic differentiation}
\item{DeltaD.per.loc}{Delta D per locus}

%% ...
}
\references{
Felsenstein, J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. Evolution, 39(4), 783-791.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
boot.ppDeltaD(infile, ncode=3, nboot=999, quant = c(0.025, 0.975), diploid = TRUE)
}

