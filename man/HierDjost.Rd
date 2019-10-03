\name{HierDjost}
\alias{HierDjost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hierarchical DJost
%%  ~~function to do ... ~~
}
\description{ This function partitions the DJost into hierarchical levels
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HierDjost(x, nreg, r, ncode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Input file, genepop format
%%     ~~Describe \code{x} here~~
}
  \item{nreg}{ The number of aggregates
%%     ~~Describe \code{nreg} here~~
}
  \item{r}{ The number of subaggregates in each aggregates
%%     ~~Describe \code{r} here~~
}
  \item{ncode}{ The integer used to code an allele
%%     ~~Describe \code{ncode} here~~
}
}
\details{ This function partitions genetic differentiation, DJost, into different hierarchical levels
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{HierDjost_perloc}{DJost per locus}
\item{HieDjost_overall}{Overall DJost}
\item{pairDjostp}{Pairwise DJost between subaggregates (populations)}
\item{pairDjostr}{Pairwise DJost between aggregates (regions)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{HierDjost_perloc}{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Jost, L. O. U. (2008). GST and its relatives do not measure differentiation. Molecular ecology, 17(18), 4015-4026.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
HierDjost_Island=HierDjost(infile,nreg=4,r=c(7,4,2,3),ncode=3)
print(HierDjost_Island)
}


