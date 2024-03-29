#' fda.oce: A package for the use of functional data analysis tools on oceanographic profiles.
#'
#' The \pkg{fda.oce} package provides a serie of functions companion to Pauthenet et al. (2017), to analyze the variability in multivariate hydrographic profiles.
#'
#' @section Functions :
#' \code{bspl} : function to fit B-splines on multivariate hydrographic profiles
#'
#' \code{fpca} : compute a basis of modes from multivariate hydrographic profiles using a functional principal component analysis (fpca).
#'
#' \code{proj} : project multivariate hydrographic profiles on a basis (computed with \code{\link{fpca}}) to obtain the principal components (pc) of the dataset
#'
#' \code{reco} : reconstruct the profiles with a chosen numbers of principal components
#'
#' @section Visualisations :
#' \code{eigenf} : effect of each eigenfunction (i. e. vertical mode) on the mean profile.
#'
#' \code{pc_plot} : plot the PC values with the appropriate axis and percentage values on the labels
#'
#'
#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#' @docType package
#' @name fda.oce
NULL
