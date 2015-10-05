\encoding{UTF-8}
\name{Dirichlet}
\alias{rdirichlet}
\title{The Dirichlet Distribution}
\description{Density function and random number generation for the Dirichlet distribution}
\usage{
rdirichlet(n, alpha)
}
\arguments{
  \item{n}{number of random observations to draw}
  \item{alpha}{the Dirichlet distribution's parameters. Can be a vector (one set of parameters for all observations) or a matrix (a different set of parameters for each observation), see \dQuote{Details}}
}
\details{
  The Dirichlet distribution is a multidimensional generalization of the Beta distribution where each dimension is governed by an \eqn{\alpha}{alpha}-parameter.
  Formally this is
  \deqn{%
    \mathcal{D}(\alpha_i)=\left[\left.\Gamma(\sum_{i}\alpha_i)\right/\prod_i\Gamma(\alpha_i)\right]\prod_{i}y_i^{\alpha_i-1}%
  }{%
    D(\alpha)=[\Gamma(\sum\alpha)/\prod\Gamma(\alpha)]\prod y^(\alpha-1)%
  }
  
  Usually, \code{alpha} is a vector thus the same parameters will be used for all observations.
  If \code{alpha} is a matrix, a complete set of \eqn{\alpha}{alpha}-parameters must be supplied for each observation.
}
\value{
returns a matrix with random numbers according to the supplied alpha vector or matrix.
}
\author{Chong Wu}
\examples{
X1 <- rdirichlet(100, c(5, 5, 10))
X1
}
