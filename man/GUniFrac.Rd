\name{GUniFrac}
\alias{GUniFrac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized UniFrac distances for comparing microbial communities.
}
\description{
A generalized version of commonly used UniFrac distances. It is defined as:
\deqn{d^{(\alpha)} = \frac{\sum_{i=1}^m b_i (p^A_{i} + p^B_{i})^\alpha 
\left\vert \frac{ p^A_{i} - p^B_{i} }{p^A_{i} + p^B_{i}} \right\vert } { 
 \sum_{i=1}^m b_i (p^A_{i} + p^B_{i})^\alpha},}{d^{(\alpha)} = 
 [\sum_{i=1}^m b_i (p^A_{i} + p^B_{i})^\alpha 
|p^A_{i} - p^B_{i}|/(p^A_{i} + p^B_{i})]/ [ 
 \sum_{i=1}^m b_i (p^A_{i} + p^B_{i})^\alpha],}
 
where \eqn{m}{m} is the number of branches, \eqn{b_i}{b_i} is the length of 
\eqn{i}{i}th branch, \eqn{p^A_{i}, p^B_{i}}{p^A_{i}, p^B_{i}} are the branch
proportion for community A and B. 

Generalized UniFrac distance contains an extra parameter \eqn{\alpha}{\alpha} 
controlling the weight on abundant lineages so the distance is not dominated 
by highly abundant lineages. \eqn{\alpha=0.5}{\alpha=0.5} has overall the best 
power. 

}
\usage{
GUniFrac(otu.tab, tree,alpha = c(0,0.5,1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{otu.tab}{
%%     ~~Describe \code{otu.tab} here~~
OTU count table, row - n sample, column - q OTU
}
  \item{tree}{
%%     ~~Describe \code{tree} here~~
Rooted phylogenetic tree of R class \dQuote{phylo}
}
  \item{alpha}{Parameter controlling weight on abundant lineages}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Return a list containing
		\item{d0}{UniFrac(0)}
		\item{d5}{UniFrac(0.5)}
		\item{d1}{UniFrac(1), weighted UniFrac}
    or a list containing
    \item{GUniFrac}{The distance matrix for different alpha}
    \item{alpha}{The weight}
}
\references{
%% ~put references to the literature/web site here ~
Chen, Jun, et al (2012). "Associating microbiome composition with environmental covariates using generalized UniFrac distances." Bioinformatics 28(16):2106-2113.
}
\author{
%%  ~~who you are~~
Chong Wu <chongwu@umn.edu>
}
\note{
%%  ~~further notes~~
The time consuming part is written in C and faster than the original one. The function only accepts rooted tree. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

groups <- throat.meta$SmokingStatus

# Calculate the UniFracs
unifracs <- GUniFrac(throat.otu.tab, throat.tree)
unifracs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distance }
\keyword{ UniFrac }
\keyword{ ecology }% __ONLY ONE__ keyword per line
