# MiSPU
Microbiome Based Sum of Powered Score (MiSPU) Tests 

## Description
**Motivation**: There is an increasing interest in investigating how the compositions of microbial communities are associated with various risk factors like environmental exposures. Distance-based analysis and micro- biome regression based kernel association test (MiRKAT) are two popular methods for such analyses. A proper choice of a phylogenetic distance is critical for the power of these methods. However, existing phylogenetic distance metrics are designed without accounting for differential information contents with various microbial lineages: given that not all microbial lineages are expected to be associated with the risk factor of interest, using all the lineages in distance calculations introduces noises, leading to power loss in the subsequent association testing.

**Results**: We propose a class of microbiome based sum of powered score (MiSPU) tests based on a newly defined generalized taxon proportion that combines observed microbial composition information with phylogenetic tree information. Different from the existing methods, a MiSPU test is based on a weighted score of the generalized taxon proportion in a general framework of regression, upweighting more likely to be associated microbial lineages. Our simulations demonstrated that one or more MiSPU tests were more powerful than MiRKAT while correctly controlling type I error rates. An adaptive MiSPU (aMiSPU) test is proposed to combine multiple MiSPU tests with various weights, approximating the most powerful MiSPU for a given scenario, consequently being highly adaptive and high powered across various scenarios. We applied MiSPU and aMiSPU to a throat microbiome dataset, showing that microbial communities were associated with the smoking status while adjusting for potential confounders.

## Install the package

We need install the dependent packages first. To download the MiSPU successfully, we must install "devtools".

```
install.packages("vegan")
install.packages("ape")
install.packages("aSPU")
install.packages("ade4")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("vegan") # install the dependent packages from CRAN

install.packages("devtools")

library(devtools)
install_github("ChongWu-Biostat/MiSPU") # install the MiSPU packages
```

## MiSPU
MiSPU performs MiSPU and aMiSPU for microbiome data set. We use a real data example here to illustrate the performance of MiSPU and aMiSPU. Note that to save time, we set permutation time equals 1000. In the paper, we set it equals 100,000.
```
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

out = MiSPU(Y,X, throat.tree,cov,resample = "perm", model =  "binomial", pow = c(2:8, Inf), n.perm = 1000)
out

proc.time() - strat.time
```







