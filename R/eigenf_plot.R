#' Plot of the vertical modes
#'
#' This function plots the effect of each vertical mode (i. e. eigenfunction) on the mean profile.
#'
#' @param pca the list produced by the function \code{fpca}.
#' @param te Which PC to plot.
#' @param sign The PCs are invariant by their sign, so the choice of the sign depends on the feature to represent. The parameter \code{sign} can contain a vector of \code{1} or \code{-1} to inverse the sign of the PCs if wanted. \code{sign} can also be used as a factor to increase the effect of eigenfunctions and see better the small variations.
#'
#' @return  Plot of the effects of the vertical mode \code{te} on the mean profiles. The curves show the mean profile (solid) and the effects of adding (red) and subtracting (blue) eigenfunctions. The percentages next to the header titles are the amount of variance explained by the mode displayed. The percentages in the horizontal axis label are the variance contained by each variable (T and S) on the mode displayed. \code{sign} can be used as a factor to enhance the deformation of + and - curves.

#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{proj}} for computing Principal Components, \code{\link{reco}} for reconstructing profiles with less modes.

#' @export
eigenf_plot <- function(pca,te,sign = 1){
  nbas  = pca$nbas
  ndim   = dim(coef)[3]
  basis = pca$basis
  Cm    = pca$Cm
  prange = basis$rangeval
  depth  = prange[1]:prange[2]

  par(mfrow=c(1,ndim), oma = c(0,0,2,0), mar=c(4.1,4.1,0.5,0.35))
  for(k in 1:ndim){
    d  = ((k-1)*nbas+1):(k*nbas)                    # Position of the variable k in the eigenvectors
    pb = round(100*sum(pca$vecnotWM[d,te]^2),0)    # Percentage of the bloc
    fdobj_vec = fda::fd(pca$axe[d,te],pca$basis,pca$fdnames)
    fdobj_m   = fda::fd(Cm[d],pca$basis,pca$fdnames)
    vp = fda::eval.fd(depth,fdobj_m + sign*fdobj_vec)
    vm = fda::eval.fd(depth,fdobj_m - sign*fdobj_vec)
    prof_m = fda::eval.fd(depth,fdobj_m)

    plot(prof_m, depth, typ = 'l'
      ,xlab = paste(pca$fdnames[[k+2]]," (",pb," %)",sep="")
      ,ylab = pca$fdnames[[1]]
      ,xlim = range(vp,vm,prof_m)
      ,ylim = c(prange[2],prange[1])
      ,yaxs = "i",las = 1)
    lines(vp, depth, col = 2)
    lines(vm, depth, col = 4)
  }
  mtext(paste("PC",te," (",round(pca$pval[te],0)," %)",sep=""),outer = TRUE,cex = 1.5)
}
