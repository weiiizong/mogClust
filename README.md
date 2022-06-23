# Multivariate Outcome-guided Clustering (mogClust)

An unified model for identifing molecular subtypes that are relevant to multiple model-selected outcomes collectively.

## Installation
* In R console

```{R}
library(devtools)
install_github("https://github.com/weiiizong/mogClust")
```

* Alternatively, download the tar.gz zipped file and install using the code below
```{R}
install.packages("~/mogClust_0.1.0.tar.gz",repos=NULL,type="source")
```

## Citation
To be updated

## Demo with LGRC lung disease data.  
* Call the LGRC lung disease data in package. The data has gene expression of top 1000 variant genes for n=259 pateints, three prognostic covariates (age, gender, BMI) and seven outomes (fev1pd1a, fvcprd1, ratiopre, RV, WBCDIFF1, WBCDIFF4, bode). See our paper for details.
* G is a gene expression matrices with 259(samples) rows and 1000(genes) columns. X is a covariate matrix with 259(samples) rows and 3(covariates) columns. Y is a matrix with 259(samples) rows and 7 continous outcomes.
* 
