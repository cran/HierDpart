\name{HierJd}
\alias{HierJd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Hierarchical Jaccard dissimilarity
%%  ~~function to do ... ~~
}
\description{ This function calculates hierarchical genetic Jaccard dissimilarity
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HierJd(f, ncode, nreg, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ Input file, genepop format
%%     ~~Describe \code{f} here~~
}
  \item{ncode}{ The integer used to code an allele
%%     ~~Describe \code{ncode} here~~
}
  \item{nreg}{ The number of aggregates (regions)
%%     ~~Describe \code{nreg} here~~
}
  \item{r}{ The number of subaggregates in each aggregate
%%     ~~Describe \code{r} here~~
}
}
\details{ This function partitions genetic Jaccard dissimilarity into different hierarchical levels
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Jdpop }{Jaccard dissimilarity between subaggregates}
    \item{Jdr }{Jaccard dissimilarity between aggregates}
    \item{HierJd }{The partitioning of Jaccard dissimilarity at different hierarchical levels}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Jaccard, P. (1901). Etude comparative de la distribution florale dans une portion des Alpes et des Jura. Bull Soc Vaudoise Sci Nat, 37, 547-579.
Jaccard, P. (1912). The distribution of the flora in the alpine zone. 1. New phytologist, 11(2), 37-50.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
HierJd_Island=HierJd(infile,nreg=4,r=c(7,4,2,3),ncode=3)
print(HierJd_Island)

}

