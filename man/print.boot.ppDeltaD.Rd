\name{print_boot_ppDeltaD}
\alias{print_boot_ppDeltaD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Print the bootstrapping result for pairwise allelic differentiation (Delta D)
%%  ~~function to do ... ~~
}
\description{Print the confidence interval of pairwise allelic differentiation (Delta D)
}
\usage{
print_boot_ppDeltaD(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ boot.ppDeltaD object
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{ Print the bootstrapping result for pairwise allelic differentiation (Delta D)
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{value}{print the confidence interval of the pairwise allelic differentiation}
}

\examples{
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
boottest=boot.ppDeltaD(infile, ncode=3, nboot=999, quant = c(0.025, 0.975), diploid = TRUE)
print_boot_ppDeltaD(boottest)
}

