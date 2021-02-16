#' Functional Principal Component Analysis of multivariate hydrographic profiles
#'
#' This function computes a FPCA on multivariate hydrographic profiles.
#' It returns the basis of decomposition, also called the vertical modes.
#'
#' @param fdobj fd objects (list) of the splines construction containing coefficients, etc... This is produced by the function \code{bspl}.

#' @return \code{pca} list containing the eigen (\code{values}), eigen (\code{vectors}), eigen vectors not weighted by WM (\code{vecnotWM}), deformation induced by an axis (\code{axe}) percentage of each axis (\code{pval}).
#'
#' @author Etienne Pauthenet \email{<etienne.pauthenet@gmail.com>}, David Nerini \code{<david.nerini@univ-amu.fr>}, Fabien Roquet \code{<fabien.roquet@gu.se>}, Alejandro Ariza

#' @references Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, https://doi.org/10.1175/JPO-D-16-0083.1
#' @references Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
#'
#' @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{proj}} for computing Principal Components, \code{\link{reco}} for reconstructing profiles with less modes.


#' @export
fpca <-function(fdobj){
  coef  = fdobj$coefs
  basis = fdobj$basis
  metric = fda::eval.penalty(basis)
  nbas   = basis$nbasis
  nobs   = dim(coef)[2]
  ndim   = dim(coef)[3]
  prange = basis$rangeval
  depth  = prange[1]:prange[2]

  C = NULL
  for(k in 1:ndim){
    C  <- cbind(C,t(coef[,,k]))
  }
  Cm <- apply(C,2,mean,na.rm = T)
  Cc <- sweep(C,2,Cm,"-")

  #Inertia
  iner = NULL
  V        <- 1/nobs*t(Cc[,1:nbas])%*%Cc[,1:nbas]%*%metric
  iner[1]  <- sum(diag(V))
  for(n in 2:ndim){
    V        <- 1/nobs*t(Cc[,(nbas*(n-1)+1):(nbas*(n-1)+nbas)])%*%Cc[,(nbas*(n-1)+1):(nbas*(n-1)+nbas)]%*%metric
    iner[n]  <- sum(diag(V))
  }

  #Metric W
  W <- matrix(0, nbas*ndim, nbas*ndim)
  for (i in 1:ndim){
    i0              <- i*nbas - nbas + 1
    i1              <- i*nbas
    W[i0:i1, i0:i1] <- metric
  }
  W       <- (W + t(W))/2
  Wdem    <- chol(W)
  Wdeminv <- solve(Wdem)

  #Metric M
  M <- rep(1/iner[1],nbas)
  for(n in 2:ndim){
    M <- c(M,rep(1/iner[n],nbas))
  }
  Mdem    <- diag(sqrt(M))
  Mdeminv <- diag(1/sqrt(M))
  M       <- Mdem^2

  ###Variance Covariance matrix weighted by WM
  V   <- 1/nobs*Mdem%*%Wdem%*%t(Cc)%*%Cc%*%t(Wdem)%*%Mdem
  pca <- eigen(V)

  #make sure that the small eigen values are not slightly negative
  pca$values <- abs(pca$values)

  ###Normalisation of eigen vectors for the metric WM
  pca$vecnotWM <- pca$vectors
  pca$vectors  <- Mdeminv%*%Wdeminv%*%pca$vectors

  ###Eigen function
  pca$values[which(Re(pca$values)<10^-10)] <- 0   #remove the zero machine
  pca$axe <- sweep(pca$vectors,2,sqrt(pca$values),"*")

  ###Percentage of each axis
  pca$pval <- Re(round(pca$values/sum(pca$values)*100,2))

  ###Store the inertia and the mean profile in the list pca, to be used for projection (proj)
  pca$iner <- iner
  pca$Cm    <- Cm

  pca$basis = basis
  pca$metric = metric
  pca$W = W
  pca$M = M
  pca$nbas = nbas
  pca$nobs = nobs
  pca$ndim = ndim
  pca$fdnames = fdobj$fdnames

  ####Verification
  # Eigen vectors should be orthogonals for the metric WM
  Verif0 <- as.numeric(t(pca$vect[,1])%*%W%*%M%*%pca$vect[,2])      # This should be = 0 (<10^-14)
  Verif1 <- as.numeric(t(pca$vect[,1])%*%W%*%M%*%pca$vect[,1] - 1)  # This should be = 0  (<10^-14)
  if(Verif0>10^-14 | Verif1 >10^-14){
    cat("Warning : the modes are not orthogonals.")
  }
  #Return
  pca <<- pca
}
