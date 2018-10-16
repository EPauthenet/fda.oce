# Functional Data Analysis of oceanographic profiles (fda.oce)

**Functional Data Analysis** is a set of tools to study curves or functions. Here we see vertical profiles of temperature (T) and salinity (S) as curves and apply a functional principal component analysis (FPCA) in the multivaraite case (T and S) to reduce the dimensionality of the system. It can be used for front detection, water mass identification, unsupervised or supervised classification, model comparison, data calibration ...

*Reference*: 
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
Here is an example of how to use these functions. We compute the modes for a subsample (given here in a RData) of the reanalysis [GLORYS](http://marine.copernicus.eu/services-portfolio/access-to-products/) in the Southern Ocean for December of 2015.

First we load the data and fit the Bsplines on the 1691 profiles of the example, between 5 and 1000m with 20 Bsplines :
``` r
load("GLORYS_2015-12_SO_sub10.RData")
fda.oce::bspl(temp,sal,depth,range = c(5,1000),mybn = 20)
```

Then we apply the PCA on these Bsplines :
``` r
fda.oce::fpca(temp.fd,sal.fd)
```

The profiles can then be reconstructed with less PCs than the total number, removing the samll variability :
``` r
fda.oce::reco(pca,te)
```

And we can project any other profile on the modes computed with a climatology (here a subset of GLORYS 12-2015). For example we fitted Bsplines on the [WOCE section i06s](https://www.nodc.noaa.gov/woce/wdiu/) and store it in a RData. Now we load it and project it on the GLORYS modes :

``` r
load("i06s_1000m40BS.RData")
fda.oce::proj(ti06.fd,si06.fd,pca)
```


# Visual representation
We can plot the effects of the vertical modes on temperature and salinity mean profiles :
``` r
for (te in 1:4){
  fda.oce::eigenf(pca,te)
}
```


<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen1.png" alt="drawing" width="350px"/> <img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen2.png" alt="drawing" width="350px"/> <img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen3.png" alt="drawing" width="350px"/> <img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen4.png" alt="drawing" width="350px"/>


We can plot the PC in space, here are the 4 first :
``` r
par(mfrow = c(2,2),mar = c(1,3,3,3))
for (te in 1:4){
  fda.oce::map_pc(lon,lat,pca,mask,te)
}
```
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_PCmap.png" alt="drawing" width="1000px"/>

We can also make a kernel density estimation (KDE) of the PC, here PC1 against PC2, with the section WOCE i06s in red :

``` r
fda.oce::kde_pc(pca, te = c(1, 2))
points(Npc[,1:2],col = 2,pch = "+")
legend("topright",legend = "section WOCE i06s",pch = "+",col = 2)
```

<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_pca.png" alt="drawing" width="1000px"/>

# Data Validation
Here is an example of data calibration on a subset of [ARGO floats](http://www.seanoe.org/data/00311/42182/) in the Southern Ocean. Any profiles that is not located in the main cloud of points can be easily flagged as of low quality.

<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/Argo_outliers.png" alt="drawing" width="1000px"/>

