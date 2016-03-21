\name{MiSPU-package}
\alias{MiSPU-package}
\docType{package}
\title{
Microbiome Based Sum of Powered Score (MiSPU) Tests 
}
\description{
There is an increasing interest in investigating how the compositions of microbial communities are associated with human health and disease. In this package, we present a novel global testing method called aMiSPU, that is highly adaptive and thus high powered across various scenarios, alleviating the issue with the choice of a phylogenetic distance. Our simulations and real data analysis demonstrated that aMiSPU test was often more powerful than several competing methods while correctly controlling type I error rates.
}
\details{
\tabular{ll}{
Package: \tab MiSPU\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-02-29\cr
License: \tab GPL-2 \cr
}
}
\author{
Chong Wu, Wei Pan
Maintainer: Chong Wu <wuxx0845@umn.edu>
}
\references{
Pan, W., et al.(2014) A powerful and adaptive association test for rare variants, Genetics, 197(4), 1081-95

Chong, W., Pan, W. (2015) An Adaptive Association Test for Microbiome Data, submitted. 
}

\keyword{MiSPU-package}

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

out = MiSPU(Y,X, throat.tree,cov,model =  "binomial", pow = c(2:8, Inf), n.perm = 1000)
out
}
