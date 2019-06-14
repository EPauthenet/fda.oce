#' Graphical representation of the PC on the main principal plans
#'
#' This function plots the PC values with the appropriate axis and percentage values on the labels.
#'
#' @param \code{pca} list containing the vertical modes used (computed with the function \code{fpca}).
#' @param \code{pc} The principal components of the profiles (computed with the function \code{proj}).
#' @param \code{n} vector of the PC to plot, it can be of length 2 or 3. Default is n = c(1,2)

#' @return Plot the main PC with their corresponding variance on the axis labels.

#'
#'#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{proj}} for computing Principal Components, \code{\link{reco}} for reconstructing profiles with less modes.


#' @export
pc_plot <-function(pca,pc,n = c(1,2)){
  if(length(n)==2){
    plot(pc[,n],las = 1,pch = 20,cex = .4
      ,xlab=paste("PC",n[1]," (",pca$pval[n[1]]," %)",sep = "")
      ,ylab=paste("PC",n[2]," (",pca$pval[n[2]]," %)",sep = ""))
    abline(v=0,h=0,lty=3)
  }
  if(length(n)==3){
    rgl::plot3d(pc[,n],col = 1,pch = 20
    ,xlab=paste("PC",n[1]," (",pca$pval[n[1]]," %)",sep = "")
    ,ylab=paste("PC",n[2]," (",pca$pval[n[2]]," %)",sep = "")
    ,zlab=paste("PC",n[3]," (",pca$pval[n[3]]," %)",sep = ""))
  }
}




