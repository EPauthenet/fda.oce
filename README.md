# Functional Data Analysis of oceanographic profiles (fda.oce)

**Functional Data Analysis** is a set of tools to study curves or functions. Here we see vertical hydrographic profiles of several variables (temperature, salinity, oxygen,...) as curves and apply a functional principal component analysis (FPCA) in the multivaraite case to reduce the dimensionality of the system. The classical case is done with couples of temperature and salinity. It can be used for front detection, water mass identification, unsupervised or supervised classification, model comparison, data calibration ...

*References*: 
- Pauthenet et al. (2018) Seasonal meandering of the Polar Front upstream of the Kerguelen Plateau. Geophysical Research Letters, [10.1029/2018GL079614](https://doi.org/10.1029/2018GL079614)
- Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, [10.1175/JPO-D-16-0083.1](http://dx.doi.org/10.1175/JPO-D-16-0083.1)
- Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.

Installation of the package using devtools :
``` r
install.packages("devtools")
devtools::install_github("Epauthenet/fda.oce")
help(package = fda.oce)
```

# Demo
Here is an example of how to use these functions. We compute the modes for temperature and salinity profiles of the reanalysis [GLORYS](http://marine.copernicus.eu/services-portfolio/access-to-products/) in the Southern Ocean for December of 2015.

First we load the data and fit the Bsplines on the 1691 profiles of the example. By default the function fit 20 Bsplines. It returns a fd object named 'fdobj' :
``` r
load("GLORYS_SO_2015-12.RData")
fda.oce::bspl(Pi,Xi)
```

Then we apply the FPCA on the fd object :
``` r
fda.oce::fpca(fdobj)
```

The profiles can be projected on the modes defined by the FPCA, to get the principal components (PCs) :
``` r
fda.oce::proj(fdobj,pca)
```

Visualisation of the 2 first PCs :
``` r
fda.oce::pc_plot(pca,pc,c(1,2))
```
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/pc_plot.png" alt="drawing" width="1000px"/>

Visualisation of the 2 first eigenfunctions effect on the mean profile (red (+1) and blue (-1)) :
``` r
fda.oce::eigenf_plot(pca,1)
fda.oce::eigenf_plot(pca,2)
```
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/eigenf1.png" alt="drawing" width="370px"/> <img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/eigenf2.png" alt="drawing" width="370px"/>

The profiles can then be reconstructed with less PCs than the total number, removing the small variability. For example with only 5 modes :
``` r
te = 5
fda.oce::reco(pca,pc,te)
```

To transform fd objects back in a the variable space, we use the function eval.fd ("fda" package) :
``` r
X = eval.fd(Pi,fdobj)
X_reco = eval.fd(Pi,fdobj_reco)
```

And finally we can represent the profiles reconstructed compared to the original data :
``` r
i = 600  #index of a profile
par(mfrow = c(1,pca$ndim))
for(k in 1:ndim){ #Loop for each variable                        
  plot(Xi[,i,k],Pi,las = 1,cex = .2,col = 1   #Plot of the raw data
    ,xlim = range(Xi[,i,k],data_reco[,i,k])
    ,ylim = c(1000,0)
    ,xlab = pca$fdnames[[2+k]]
    ,ylab = pca$fdnames[[1]])
  lines(X[,i,k],Pi,col = 2)              #Plot of the B-spline fit
  lines(X_reco[,i,k],Pi,col = 3)         #Plot of the reconstructed profiles
}
legend("bottomleft",col = c(1,2,3),lty = c(NA,1,1),pch = c(20,NA,NA)
  ,legend = c("raw data","B-spline fit",paste("reconstruction with ",te," modes",sep = "")))
```
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/reco_prof600.png" alt="drawing" width="1000px"/>

