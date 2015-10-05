\name{MiSPU-package}
\alias{MiSPU-package}
\alias{MiSPU}
\docType{package}
\title{
\packageTitle{MiSPU}
}
\description{
\packageDescription{MiSPU}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MiSPU}
\packageIndices{MiSPU}}
\author{
\packageAuthor{MiSPU}

Maintainer: \packageMaintainer{MiSPU}
}
\references{
Pan, W., et al.(2014) A powerful and adaptive association test for rare variants, Genetics, 197(4), 1081-95

Chong, W., Pan, W. (2015) An Adaptive Association Test for Microbiome Data, submitted. 
}

\keyword{ package }

\examples{
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

Y.tmp =throat.meta[,3]
Y = rep(0,dim(throat.meta)[1])
Y[Y.tmp=="Smoker"] = 1

cov.tmp = throat.meta[,c(10,12)]
cov = matrix(1,dim(throat.meta)[1],2)
cov[cov.tmp[,1]== "None",1] = 0
cov[cov.tmp[,2]== "Male",2] = 0

start.time = proc.time()
X = as.matrix(throat.otu.tab)

out = MiSPU(Y,X, throat.tree,cov,resample = "perm",
      model =  "binomial", pow = c(2:8, Inf), n.perm = 1000)
out
}
