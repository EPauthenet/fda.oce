#' Projection of T-S profiles
#'
#' Projection of any temperature and salinity (T-S) profiles on a chosen basis. It is essential to locate a profile relatively to a climatology. The profiles to project must have the same \code{myb} than the modes to project on (i.e. same length of profile \code{c(dmin,dmax)} and same number of Bsplines \code{mybn}).
#'

#' @param temp.fd,sal.fd fd objects (list) of the T-S profile(s) to project. This is produced by the function \code{bspl}.
#' @param pca list of the modes to project on, produced by the function \code{fpca}.
#' @param te how many PCs to use in the reconstruction, default is set to the total number of PC, \code{te = mybn}.

#' @return \code{Npc} The principal components of the profiles \code{temp.fd} and \code{sal.fd} projected on the modes contained in \code{pca}.
#'
#'@references Ramsay J. O., and B. W. Silverman, 2005: Functional Data Analysis. Springer, 426 pp.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} FPCA of T-S profiles, \code{\link{PCmap}} for plotting a map of PC, \code{\link{kde_pc}} for kernel density estimation of two PCs...

proj <- function(temp.fd,sal.fd,pca){
  mybn  <- pca$nbasis
  myb   <- create.bspline.basis(c(pca$range[1],pca$range[2]),mybn)
  Cm = pca$Cm

  #Metric M
  M=c(rep(1/pca$inerT,mybn),rep(1/pca$inerS,mybn))
  Mdem=diag(sqrt(M))
  Mdeminv=diag(1/sqrt(M))
  M=Mdem^2

  #Metric W
  W=eval.penalty(myb)
  nul=matrix(0,mybn,mybn)
  W=cbind(rbind(W,nul),rbind(nul,W))
  W=(W+t(W))/2

  #Combine T and S and substract the mean of the modes to project on
  C=cbind(t(temp.fd$coefs),t(sal.fd$coefs))
  Cc=sweep(C,2,Cm,"-")

  Npc<<-Cc%*%t(W)%*%M%*%pca$vectors
}
