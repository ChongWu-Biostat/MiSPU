# MiSPU
Microbiome Based Sum of Powered Score (MiSPU) Tests 

## Install the package
We test it on R 3.2.1 in Linux server and 3.2.2 in Windows and Mac. *Note that for windows user, we need install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. For Mac user, we need install [gfortran](https://cran.r-project.org/bin/macosx/tools/).*
```
install.packages("vegan") # 
install.packages("ape")
install.packages("aSPU")
install.packages("ade4")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("vegan") # install the dependent packages from CRAN.

install.packages("devtools")

library(devtools)
install_github("ChongWu-Biostat/MiSPU") # install the MiSPU packages
```

## MiSPU
MiSPU performs MiSPU and aMiSPU for microbiome data set. We use a real data example here to illustrate the performance of MiSPU and aMiSPU. Note that to save time, we set permutation time equals 1000. In the paper, we set it equals 100,000.
```
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

out = MiSPU(Y,X, throat.tree,cov,resample = "perm", model =  "binomial", pow = c(2:8, Inf), n.perm = 1000)
out

proc.time() - start.time
```

## Generalized UniFrac
We note that calculating generalized unifrac distance matrix is time consuming especially for the large data set. To save computational time, here, we provide a C version of generalized UniFrac function. We expect this function is much faster than the existing GUniFrac R package.
```
data(throat.otu.tab)
data(throat.tree)
data(throat.meta) # the data is from the paper: Charlson ES, Chen J, Custers-Allen R, Bittinger K, Li H, et al. (2010) Disordered Microbial Com- munities in the Upper Respiratory Tract of Cigarette Smokers. PLoS ONE 5(12): e15216.


groups <- throat.meta$SmokingStatus

# Calculate the UniFracs
unifracs <- GUniFrac(throat.otu.tab, throat.tree)
unifracs
```

## Manual
If you like MiSPU, please give us a star. You can download the MiSPU package manual [here](http://cutpi.com/upimages/1444077574.pdf). 






