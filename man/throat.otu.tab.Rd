\name{throat.otu.tab}
\alias{throat.otu.tab}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
OTU count table from 16S sequencing of the throat microbiome samples.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
It is part of a microbiome data set for studying the effect of
 smoking on  the upper respiratory tract microbiome. The original data 
set contains samples from both throat and  nose microbiomes, and 
from both body sides. This data set comes from the throat microbiome 
of left body side. It contains 60 subjects consisting of 32 nonsmokers
 and 28 smokers. 
}
\usage{data(throat.otu.tab)}
\format{
  The format is:
 chr "throat.otu.tab"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The OTU table is produced by the QIIME software. Singleton OTUs have
been discarded.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Charlson ES, Chen J, Custers-Allen R, Bittinger K, Li H, et al. (2010) 
Disordered Microbial Communities in the Upper Respiratory Tract of 
Cigarette Smokers. PLoS ONE 5(12): e15216. 
}

\examples{
data(throat.otu.tab)
## maybe str(throat.otu.tab) ; plot(throat.otu.tab) ...
}
\keyword{datasets}
