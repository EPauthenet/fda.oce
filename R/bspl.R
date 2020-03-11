#' B-spline fits on Multivariate Hydrographic Profiles
#'
#' This function fits B-splines on multivariate hydrographic profiles and return a functional data object.
#'
#' @param Xi array containing the profiles stored in this order \code{levels} x \code{stations} x \code{variables}
#' @param Pi vector containing the levels
#' @param nbas number of Bsplines (coefficients), by default nbas = 20.
#' @param fdn a list of the variable names, by default fdn = list('Temperature','Salinity').
#'
#' @return \code{fdobj} : fd objects of the splines construction containing coefficients, basis etc... The coefficients are stored in an array \code{nbasis} x \code{stations} x \code{variables}

#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#'
#' @seealso \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{proj}} for computing Principal Components, \code{\link{reco}} for reconstructing profiles with less modes.

#' @export
bspl <- function(Pi,Xi,nbas = 20,fdn = list('Temperature','Salinity')){
  cat("Converting the profiles into Bsplines...")
  require(fda)
  fdn1 = list('Level','Station')
  fdnames = c(fdn1,fdn)
  prange = c(Pi[1],Pi[length(Pi)])
  Breaks=prange[1]+(prange[2]-prange[1])*tan(seq(0,1,1/(nbas-3)))/tan(1)
  basis = fda::create.bspline.basis(rangeval = prange,nbasis = nbas,norder = 4,breaks = Breaks)

  fdobj <<- fda::Data2fd(argvals = Pi,y = Xi,basisobj = basis,fdnames = fdnames)

  cat(paste(dim(Xi)[2],"B-splines computed for",dim(Xi)[3],"variables."))
}

