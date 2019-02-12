#' Functional Principal Component Analysis of Temperature and Salinity profiles
#'
#' Functional Principal Component Analysis (FPCA) in the multivariate case, applied on Temperature and Salinity (T-S) profiles seen as curves (Bsplines)
#'
#' @param temp.fd,sal.fd fd objects (list) of the splines construction containing coefficients, etc... This is produced by the function \code{bspl}.
#' @param plot,plot3d if TRUE, plot the first two or three PCs.
#' @param we to do a weighted PCA. Vector of length N containing the weights to assign to each profile.

#' @return \code{pca} list containing the eigen (\code{values}), eigen (\code{vectors}), eigen vectors not weighted by WM (\code{vecnotWM}), principal components (\code{pc}), deformation induced by an axis (\code{axe}) percentage of each axis (\code{pval}).
#'
#'@references Ramsay J. O., and B. W. Silverman, 2005: Functional Data Analysis. Springer, 426 pp.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{PCmap}} for plotting a map of PC, \code{\link{kde_pc}} for kernel density estimation of two PCs...

#' @export
fpca <-function(temp.fd,sal.fd,plot = FALSE,plot3d = FALSE,we){
  mybn  <- length(temp.fd$coefs[,1])
  dmin  <- temp.fd$basis$rangeval[1]
  dmax  <- temp.fd$basis$rangeval[2]
  nobs  <- dim(temp.fd$coefs)[2]
  myb   <- fda::create.bspline.basis(c(dmin,dmax),mybn)
  depth <- dmin:dmax

  C  <- cbind(t(temp.fd$coefs),t(sal.fd$coefs))
  Cm <- apply(C,2,mean)
  Cc <- sweep(C,2,Cm,"-")

  #Inertia
  metric <- fda::eval.penalty(myb)
  VT     <- 1/nobs*t(Cc[,1:mybn])%*%Cc[,1:mybn]%*%metric
  inerT  <- sum(diag(VT))
  VS     <- 1/nobs*t(Cc[,(mybn+1):(2*mybn)])%*%Cc[,(mybn+1):(2*mybn)]%*%metric
  inerS  <- sum(diag(VS))

  #Metric W
  W       <- fda::eval.penalty(myb)
  nul     <- matrix(0,mybn,mybn)
  W       <- cbind(rbind(W,nul),rbind(nul,W))
  W       <- (W+t(W))/2
  Wdem    <- chol(W)
  Wdeminv <- solve(Wdem)

  #Metric M
  M       <- c(rep(1/inerT,mybn),rep(1/inerS,mybn))
  Mdem    <- diag(sqrt(M))
  Mdeminv <- diag(1/sqrt(M))
  M       <- Mdem^2

  ###Variance Covariance matrix weighted by WM
  if (exists(we)==TRUE){
    V   <- 1/nobs*Mdem%*%Wdem%*%t(we*Cc)%*%Cc%*%t(Wdem)%*%Mdem
  }else{
    V   <- 1/nobs*Mdem%*%Wdem%*%t(Cc)%*%Cc%*%t(Wdem)%*%Mdem
  }
  pca <- eigen(V)

  #make sure that the small eigen values are not slightly negative
  pca$values <- abs(pca$values)

  ##Principal Components (PCs)
  pca$pc <- Cc%*%t(Wdem)%*%Mdem%*%pca$vectors

  ###Normalisation of eigen vectors for the metric WM
  pca$vecnotWM <- pca$vectors
  pca$vectors  <- Mdeminv%*%Wdeminv%*%pca$vectors

  ###Eigen function
  pca$values[which(Re(pca$values)<10^-10)] <- 0   #remove the zero machine
  pca$axe <- sweep(pca$vectors,2,sqrt(pca$values),"*")

  ###Percentage of each axis
  pca$pval <- Re(round(pca$values/sum(pca$values)*100,2))

  ###Store the inertia and the mean profile in the list pca, to be used for projection (proj)
  pca$inerT <- inerT
  pca$inerS <- inerS
  pca$Cm    <- Cm

  ###Store the informations of the basis myb
  pca <-c(pca,myb[2:5])

  ####Verification
  # Eigen vectors should be orthogonals for the metric WM
  Verif0 <- as.numeric(t(pca$vect[,1])%*%W%*%M%*%pca$vect[,2])      # This should be = 0 (<10^-14)
  Verif1 <- as.numeric(t(pca$vect[,1])%*%W%*%M%*%pca$vect[,1] - 1)  # This should be = 0  (<10^-14)
  if(Verif0>10^-14 | Verif1 >10^-14){
    cat("Warning : the modes are not orthogonals.")
  }
  Verif3 <- pca$val[1] - 1/nobs*t(pca$pc[,1])%*%pca$pc[,1]
  if(Verif3>10^-10){
    cat("Warning : The eigen values are not equivalent to the PC norm.")
  }

  if (plot == TRUE){
    plot(pca$pc[,1:2],xlab=paste("PC1 (",pca$pval[1]," %)"),las = 1
      ,ylab=paste("PC2 (",pca$pval[2]," %)"),col = 1,pch = 20,cex = .4)
    abline(v=0,h=0,lty=3)
  }
  if (plot3d == TRUE){
    rgl::plot3d(pca$pc[,1:3],xlab=paste("PC1 (",pca$pval[1]," %)")
      ,ylab=paste("PC2 (",pca$pval[2]," %)")
      ,zlab=paste("PC3 (",pca$pval[3]," %)"),col = 1,pch = 20)
  }

  #Return
  pca <<- pca
}
