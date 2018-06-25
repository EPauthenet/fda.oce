# Functional Data Analysis of oceanographic profiles (fda.oce)

**Functional Data Analysis** is an emerging set of tools to study curves or functions. Here we see vertical profiles of temperature (T) and salinity (S) as curves and apply a functional principal component analysis (FPCA) in the multivaraite case (T and S) to reduce the dimensionality of the system. It can be used for front detection, water mass identification, unsupervised or supervised classification, model comparison, ...

*Reference*: 
- Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1


# Example
Here is an example of how to use these functions. We compute the modes for a subsample (given here in a RData) of GLORYS in the Southern Ocean for December of 2015 (mercatorglorys2v4_gl4_mean_201512.nc). The total data is available here http://marine.copernicus.eu/services-portfolio/access-to-products/ .

First we load the data and fit the Bsplines on the 1691 profiles of the example :
``` r
load("GLORYS_2015-12_SO_sub10.RData")
bspl(temp,sal,depth,range = c(5,1000),mybn = 20)
```

Then we apply the PCA on these Bsplines :
``` r
fpca(temp.fd,sal.fd,plot = T)
```

The profiles can then be reconstructed with less PCs than the total number, removing the samll variability :
``` r
reco(pca,te)
```

And we can project any other profile on the modes computed with a climatology (here a subset of GLORYS 12-2015). For example we fitted Bsplines on the WOCE section i06s (available here : https://www.nodc.noaa.gov/woce/wdiu/) and store it in a RData. Now we load it and project it on the GLORYS modes :

``` r
load("i06s_1000m40BS.RData")
proj(ti06.fd,si06.fd,pca)
```


# Visual representation
We can plot the effects of the vertical modes on temperature and salinity mean profiles :
``` r
for (te in 1:4){
  eigenf(pca,te)
}
```
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen1.png" alt="drawing" width="500px"/>
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen2.png" alt="drawing" width="500px"/>
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen3.png" alt="drawing" width="500px"/>
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_eigen4.png" alt="drawing" width="500px"/>


We can plot the PC in space, here are the 4 first :
``` r
par(mfrow = c(2,2),mar = c(1,3,3,3))
for (te in 1:4){
  map_pc(lon,lat,pca,mask,te)
}
```
<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_PCmap.png" alt="drawing" width="1000px"/>

We can also make a kernel density estimation (KDE) of the PC, here PC1 against PC2, with the section WOCE i06s in red :

``` r
kde_pc(pca, te = c(1, 2))
points(Npc[,1:2],col = 2,pch = "+")
legend("topright",legend = "section WOCE i06s",pch = "+",col = 2)
```

<img src="https://github.com/EPauthenet/fda.oce/blob/master/figures/GLO_pca.png" alt="drawing" width="1000px"/>



