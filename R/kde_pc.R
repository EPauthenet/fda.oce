#' Kernel density plot of two PCs
#'
#' Produce a kernel density estimation plot of two PCs (using \code{kde2d {MASS}}), on the axis are the percentage of explained variance of each PC (pca$pval).
#'
#' This function is made to be modified at will, to add contours, change colors etc...
#'
#' @param pca The \code{pca} list returned by the function \code{fpca}.
#' @param te A vector of which two PCs to plot against each other. The default values are PC1/PC2 (i.e. te = c(1,2)) makes a lot of sense when looking at T-S profiles.
#' @param h vector of bandwidths for x and y directions. Defaults is h = 0.2.
#' @param n Number of grid points in each direction. Can be scalar or a length-2 integer vector. Default is h = 50.
#' @param l levels of contour, default is l = seq(.1,4,.2).
#' @return A kernel density plot of two PCs i.e. a 2D boxplot.

#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{PCmap}} for plotting a map of PC

kde_pc <- function(pca,te = c(1,2),h = 0.2,n = 50,l = seq(0.1,4,.2),title = ""){
  PC1  <- pca$pc[,te[1]]
  PC2  <- pca$pc[,te[2]]
  PCm1 <- PC1[is.finite(PC1)]
  PCm2 <- PC2[is.finite(PC2)]
  d    <- kde2d(PCm1,PCm2,h = h,n = n)
  colo <- designer.colors( n=265, col= c("white","grey", "grey35","grey15"),alpha=1)

  image.plot(d,col = colo,las = 1,main = title
    ,xlab = paste('PC',te[1]," (",pca$pval[te[1]]," %)",sep = "")
    ,ylab = paste('PC',te[2]," (",pca$pval[te[2]]," %)",sep = "")
    ,xlim = range(PC1)
    ,ylim = range(PC2))
  contour(d,add = T,drawlabels = FALSE,lev = l)
  grid()
  abline(v = 0,h = 0, lty = 3)
}
