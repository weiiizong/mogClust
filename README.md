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

## Demo with LGRC lung disease data 
* Call the LGRC lung disease data in package. The data has gene expression of top 1000 variant genes for n=259 pateints, three prognostic covariates (age, gender, BMI) and seven outomes (fev1pd1a, fvcprd1, ratiopre, RV, WBCDIFF1, WBCDIFF4, bode). See our paper for details.
* G is a gene expression matrices with 259(samples) rows and 1000(genes) columns. X is a covariate matrix with 259(samples) rows and 3(covariates) columns. Y is a matrix with 259(samples) rows and 7 continous outcomes.
```{R}
X = lung_1000G$X
G = lung_1000G$G
Y = lung_1000G$Y
n = nrow(G)
np = ncol(X)
nq = ncol(Y)
NG = ncol(G)
```

* Set the number of subtypes/clusters and tuning parameters.

```{R}
K = 4
lambdaB = 0.002436381
lambdaG = 0.1507389
```

* Set initial parameters using sparse K-means.
```{R}
set.seed(12345)
library(sparcl)
fit = lm(Y~X)
C = coef(fit)[-1,]

Y1 = Y-X%*%C
km.perm = KMeansSparseCluster.permute(Y1,K=K,nperms=20)
km.out = KMeansSparseCluster(Y1,K=K,wbounds=km.perm$bestw,nstart = 150)
Cs = km.out[[1]]$Cs
w_old = t(sapply(1:length(Cs), function(x){
  a = rep(0,K)
  a[Cs[x]] = 1
  return(a)
}))
library(glmnet)
cvfit=cv.glmnet(x=G, y=w_old, family="multinomial", type.multinomial = "grouped", intercept=T, alpha = 0.5)
coef = coef(cvfit,s="lambda.1se")
gamma_new_matrix=sapply(1:K, function(x) as.numeric(coef[[x]]))[-1,]
gamma=sapply(1:K, function(x) gamma_new_matrix[,x]-gamma_new_matrix[,K])[,-K]

center.ls = lapply(1:K, function(k){
  apply(Y1[which(km.out[[1]]$Cs==k),],2,mean)
})
B = do.call(cbind,center.ls)
Sigma = diag(abs(rnorm(ncol(Y))),nrow = nq)

initials_list = list(B=B, C=C, Sigma=Sigma, gamma=gamma)

```
* Fit the mogClust model by `fit_mogClust` function.

```{R}
full_fit = fit_mogClust(G, X, Y, gamma = initials_list$gamma, Sigma = initials_list$Sigma, 
                        B = initials_list$B, C = initials_list$C, K, lambdaB = lambdaB, lambdaG.seq=lambdaG,
                        runs = 200, quite=FALSE, alphaG = 0.5, seed = "fixed")
```

* Output of the fitted model.
```{R}
str(fit,max.level = 1)
List of 12
 $ gamma        : num [1:1000, 1:3] 0 0 0 0 0 ...
 $ Sigma        : num [1:7, 1:7] 0.5657 0.4646 0.3315 -0.3312 -0.0648 ...
  ..- attr(*, "dimnames")=List of 2
 $ B            : num [1:7, 1:4] -0.368 0.74 -0.865 0.437 0 ...
 $ C            : num [1:3, 1:7] 0.8824 0.0336 0.8311 0.8079 0.0299 ...
  ..- attr(*, "dimnames")=List of 2
 $ lambdaG_path : num [1:48] 0.151 0.151 0.151 0.151 0.151 ...
 $ lambdaB      : num 0.00244
 $ R2_RMSE_table:'data.frame':	7 obs. of  2 variables:
 $ BIC          : num 4308
 $ K            : num 4
 $ pred_prob    : num [1:259, 1:4] 0.4154 0.546 0.4078 0.2951 0.0868 ...
 $ pred_lb      : Factor w/ 4 levels "1","2","3","4": 1 1 1 1 3 1 1 4 4 1 ...
 $ pred_outcome : num [1:259, 1:7] -0.6925 -0.659 -0.0102 0.0194 0.2253 ...
  ..- attr(*, "dimnames")=List of 2
```
