#' Plot of the vertical modes
#'
#' Plot the effect of each vertical mode (i. e. eigenfunction) on the mean Temperature and Salinity (T-S) profile.
#'
#'
#' @param pca the list produced by the function \code{fpca}.
#' @param te Which PC to plot.
#' @param le number of '+' or '-' by curve, default at \code{le = 50}
#' @param sign The PCs are invariant by their sign, sign can contain a vector of \code{1} or \code{-1} to inverse the sign of the PCs if wanted. \code{sign} can also be used as a factor to increase the effect of eigenfunctions and see better the small variations.
#' @param tlim,slim limits of the plot in temperature and salinity
#'
#'
#' @return  Representation of the effects of the vertical mode \code{te} on temperature and salinity mean profiles. The curves show the mean profile (solid) and the effects of adding (+) and subtracting (-) eigenfunctions. The percentages next to the header titles are the amount of variance explained by the mode displayed. The percentages in the horizontal axis label are the variance contained by each variable (T and S) on the mode displayed. \code{sign} can be used as a factor to enhance the deformation of + and - curves.

#'@references Ramsay J. O., and B. W. Silverman, 2005: Functional Data Analysis. Springer, 426 pp.

#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{PCmap}} for plotting a map of PC

#' @export
eigenf <- function(pca,te,le = 50,sign = 1,tlim = c(-2,25),slim = c(33,38)){
  if(length(sign)>1){
    s <- sign[te]
  }
  else{
    s <- 1
  }
  mybn <- pca$nbasis
  dmin <- pca$range[1]
  dmax <- pca$range[2]
  myb  <- fda::create.bspline.basis(c(dmin,dmax),mybn)
  Cm   <- pca$Cm

  #Percentage of the bloc
  pb    <- round(100*sum(pca$vecnotWM[1:mybn,te]^2),0)
  profx <- seq(dmin,dmax,length=le)

  #initialize two fd objects
  a.fd <- fda::Data2fd(dmin:dmax,rep(1,length(dmin:dmax)),myb)
  t.fd <- a.fd
  s.fd <- a.fd

  # Temperature [1:mybn,1]
  a.fd$coefs = pca$axe[1:mybn,te]
  t.fd$coefs = Cm[1:mybn]
  vpl        = (t.fd+s*a.fd)
  vmo        = (t.fd-s*a.fd)
  vplus      = fda::eval.fd(seq(dmin,dmax,length=le),vpl)
  vmoins     = fda::eval.fd(seq(dmin,dmax,length=le),vmo)
  Tm         = fda::eval.fd(seq(dmin,dmax,length=le),t.fd)

  par(mfrow=c(1,2), oma = c(0,0,2,0), mar=c(4.1,4.1,0.5,0.35))
  plot(Tm,profx,typ = 'l'
    ,xlab = paste("Temperature (C) (",pb," %)",sep=""),ylab="Depth (m)"
    ,xlim = tlim,ylim = c(dmax,dmin),yaxs = "i",las = 1)
  points(vplus,profx,pch="+")
  points(vmoins,profx,pch="-")

  # Salinity [mybn+1:2*mybn,1]
  pb    = round(100*sum(pca$vecnotWM[(mybn+1):(2*mybn),te]^2),2)
  a.fd$coefs = pca$axe[(mybn+1):(2*mybn),te]
  s.fd$coefs = Cm[(mybn+1):(2*mybn)]
  vpl        = (s.fd+s*a.fd)
  vmo        = (s.fd-s*a.fd)
  vplus      = fda::eval.fd(seq(dmin,dmax,length=le),vpl)
  vmoins     = fda::eval.fd(seq(dmin,dmax,length=le),vmo)
  Sm         = fda::eval.fd(seq(dmin,dmax,length=le),s.fd)

  par(mar=c(4.1,0.35,0.5,4.1))
  plot(Sm,profx,typ = 'l',axes = F
    ,xlab=paste("Salinity (PSU) (",pb," %)",sep="")
    ,xlim = slim,ylim = c(dmax,dmin),yaxs = "i",las = 1)
  points(vplus,profx,pch = "+")
  points(vmoins,profx,pch="-")
  axis(1)
  axis(2,labels = F)
  box()
  mtext(paste("PC",te," (",round(pca$pval[te],0)," %)",sep=""),outer = TRUE,cex = 1.5)
}
