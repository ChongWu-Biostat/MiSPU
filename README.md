# MiSPU [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![CRAN](http://www.r-pkg.org/badges/version/MiSPU)](http://cran.rstudio.com/package=MiSPU) [![Downloads](http://cranlogs.r-pkg.org/badges/MiSPU?color=brightgreen)](http://www.r-pkg.org/pkg/MiSPU)![downloads](http://cranlogs.r-pkg.org/badges/grand-total/MiSPU)
Microbiome Based Sum of Powered Score (MiSPU) Tests. The packag is avaiable at CRAN now. 

## Installation
To install the stable version from CRAN, simply run the following from an R console:

```r
install.packages("MiSPU")
```

To install the latest development builds directly from GitHub, run this instead:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("ChongWu-Biostat/MiSPU")
```

## MiSPU
MiSPU performs MiSPU and aMiSPU for microbiome data set. We use a real data example here to illustrate the performance of MiSPU and aMiSPU. Note that the number of permutation (n.perm) equals 10,000. In the paper, we set it equals 2,000,000 (Take around one or two hours to run depending on the performance of computers).
```r
library(MiSPU)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta) # the data is from the paper: Charlson ES, Chen J, Custers-Allen R, Bittinger K, Li H, et al. (2010) Disordered Microbial Com- munities in the Upper Respiratory Tract of Cigarette Smokers. PLoS ONE 5(12): e15216.

Y.tmp =throat.meta[,3]
Y = rep(0,dim(throat.meta)[1])
Y[Y.tmp=="Smoker"] = 1
cov.tmp = throat.meta[,c(10,12)]
cov = matrix(1,dim(throat.meta)[1],2)
cov[cov.tmp[,1]== "None",1] = 0
cov[cov.tmp[,2]== "Male",2] = 0
start.time = proc.time()
X = as.matrix(throat.otu.tab)

out = MiSPU(Y,X, throat.tree,cov, model =  "binomial", pow = c(2:8, Inf), n.perm = 1000)
out

proc.time() - start.time
```

## Generalized UniFrac
We note that calculating generalized unifrac distance matrix is time consuming especially for the large data set. To save computational time, here, we provide a C version of generalized UniFrac function. We expect this function is much faster than the existing GUniFrac R package. Some users stated that this function is around 4 times faster than the competing ones and we will do a fully test in the near future. 
```r
data(throat.otu.tab)
data(throat.tree)
data(throat.meta) # the data is from the paper: Charlson ES, Chen J, Custers-Allen R, Bittinger K, Li H, et al. (2010) Disordered Microbial Com- munities in the Upper Respiratory Tract of Cigarette Smokers. PLoS ONE 5(12): e15216.


groups <- throat.meta$SmokingStatus

# Calculate the UniFracs
unifracs <- GUniFrac(throat.otu.tab, throat.tree)
unifracs
```

## Manual
If you like *MiSPU*, please give us a star. You can download the *MiSPU* package manual [here](https://cran.r-project.org/web/packages/MiSPU/MiSPU.pdf). 






