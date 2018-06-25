#' fda.oce: A package for the use of functional data analysis tools on oceanographic profiles.
#'
#' The \pkg{fda.oce} package provides a serie of functions companion to Pauthenet et al. (2017), to analyze the variability in a field of temperature and salinity.
#'
#' @section Functions :
#' \code{bspl} : functions to see temperature and salinity (T-S) profiles as curves (Bsplines)
#'
#' \code{fpca} : analyze the joint variance of T and S with the use of a Functional Principal Component Analysis
#'
#' \code{reco} : The profiles can be reconstructed with less Principal Components (PCs) to see the effect of the modes separately
#'
#' \code{proj} : Any T-S profiles (observation, model,...) can be transformed into a function (\code{bspl}) and projected on the modes wanted.
#'
#' @section Visualisations :
#' \code{eigenf} : effect of each eigenfunction (i. e. vertical mode) on the mean T-S profile.
#'
#' \code{map_pc} : Map of PC
#'
#' \code{kde_pc} : Kernel density estimation of two PCs
#'
#' @author Etienne Pauthenet \code{<etienne.pauthenet@gmail.com>}, David Nerini
#'
#'
#'
#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
"_PACKAGE"
