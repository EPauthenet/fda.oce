#' Projection of hydrographic profiles
#'
#' Projection of hydrographic profiles on a chosen basis (created with \code{fpca}). It is essential to locate a profile relatively to a climatology. The profiles to project must have the same \code{basis} than the modes to project on (i.e. same length of profile \code{c(dmin,dmax)} and same number of Bsplines \code{nbasis}).
#'

#' @param fdobj fd objects (list) of the splines construction containing coefficients, etc... This is produced by the function \code{bspl}.
#' @param pca list containing the modes to project on, produced by the function \code{fpca}.

#' @return \code{pc} The principal components of the profiles stored in \code{fdobj} projected on the modes contained in \code{pca}.
#'
#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{reco}} for reconstructing profiles with less modes.

#' @export
proj <- function(fdobj,pca){
  nbas  = pca$nbasis
  basis = pca$basis
  coef  = fdobj$coefs
  Cm    = pca$Cm
  nobs  = dim(coef)[2]
  ndim   = dim(coef)[3]

  C = NULL
  for(k in 1:ndim){
    C  <- cbind(C,t(coef[,,k]))
  }
  Cc <- sweep(C,2,Cm,"-")

  pc<<-Cc %*% t(pca$W) %*% pca$M %*% pca$vectors
}



