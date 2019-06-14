#' Reconstruction of hydrographic profiles
#'
#' This function reconstructs hydrographic profiles with a chosen number of Principal Components (PCs).
#'
#' @param pca list containing the modes produced by the function \code{fpca}
#' @param Ntrunc how many PCs to use in the reconstruction, default is set to the total number of PC, \code{Ntrunc = nbas*ndim}.

#'
#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{proj}} for computing Principal Components.

#' @export
reco <- function(pca,pc,Ntrunc){
  nbas   = pca$nbas
  nobs   = dim(pc)[1]
  ndim   = dim(pc)[2]/nbas

  if(missing(Ntrunc)){Ntrunc = nbas*ndim}

  coef = array(NaN,c(nbas,nobs,ndim))
  for(k in 1:ndim){
    #coef[,,k] = repmat(pca.Cm((kk-1)*nbas+1:kk*nbas)',1,nobs) + pca.vectors((kk-1)*nbas+1:kk*nbas,1:Ntrunc)*pc(:,1:Ntrunc)'
    d = ((k-1)*nbas+1):(k*nbas)
    coef[,,k] = matrix(rep(pca$Cm[d],nobs),nbas,nobs) + pca$vectors[d,1:Ntrunc] %*% t(pc[,1:Ntrunc])
  }
    fdobj_reco <<- fd(coef,pca$basis,pca$fdnames)
}


