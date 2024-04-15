#' Map of Principal Component
#'
#' Produce a map of Principal Component (PC), with the colorbar centered in zero, the coastline, the percentage of explained variance associated with the PC (pca$pval).
#'
#' This function is made to be modified at will, to add contours, change colors etc...
#'
#' @param lon A vector of longitude.
#' @param lat A vector of latitude.
#' @param pca the list produced by the function \code{fpca}.
#' @param te Which PC is in the \code{PC} matrix.
#' @param sign The PCs are invariant by their sign, sign can contain a vector of \code{1} or \code{-1} to inverse the sign of the PCs if wanted.
#'
#' @return A color map of the PC.

#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles,...

#' @export
map_pc <- function(lon,lat,pca,mask,te,sign = 1){
  if(length(sign)>1){
    s <- sign[te]
  }

  PC                 <- mask
  PC[which(mask==1)] <- pca$pc[,te]*s

  colo <- fields::designer.colors(n=length(PC), col= c("darkblue","blue","skyblue2", "white","darkorange","red", "darkred"),alpha=1)
  zz   <- c(-max(abs(PC),na.rm = TRUE),max(abs(PC),na.rm = TRUE))

  fields::image.plot(lon,lat,PC,zlim = zz,col = colo,las = 1,xaxs = "i",yaxs = "i"
    ,main = paste('PC',te," (",pca$pval[te]," %)",sep = "")
    ,ylab = 'Longitude'
    ,xlab = 'Latitude')
  maps::map(add = T,fill = T,col = "black")
  contour(lon,lat,PC,add = TRUE)
}
