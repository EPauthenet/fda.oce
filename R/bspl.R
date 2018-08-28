#' Fit of B-splines on Temperature and Salinity profiles
#'
#' Fit B-splines on Temperature and Salinity (T-S) profiles from gridded dataset (reanalysis, models,...)

#' @param temp temperature stored in arrays \code{ind*depth} or \code{lon*lat*depth} or \code{lon*lat*depth*time}
#' @param sal salinity stored in arrays \code{lon*lat*depth} or \code{lon*lat*depth*time}
#' @param depth depth can be a vector or a matrix (case of irregular vertical sampling)
#' @param range vector of two elements containing the first and last depth level \code{c(dmin,dmax)}, \code{dmin} and \code{dmax} must be contained by the \code{depth} range.
#' @param mybn number of splines, Default is mybn = 20.
#' @param l number of interpolation point. The number of points to interpolate must be large enough to avoid overfitting. Default is l = 100.

#'
#' @return \code{temp.fd} and \code{sal.fd} : fd objects of the splines construction containing coefficients, etc...
#' @return  \code{BS} list containing :\code{profx}, \code{Teval} and \code{Seval} (the splines sampled along profx)
#' @return  \code{mask} : if the data was on a regular horizontal grid

#'
#' @seealso \code{\link{fpca}} for functional principal component analysis of T-S profiles, \code{\link{PCmap}} for plotting a map of PC, \code{\link{kde_pc}} for kernel density estimation of two PCs...

#' @export
bspl <- function(temp,sal,Depth,range,mybn = 20,l = 100){
  cat("Converting the profiles into Bsplines...")
  dmin <- range[1]
  dmax <- range[2]
  profx <- seq(dmin,dmax,length = l)

  if (length(dim(temp))==2){
    Tm <- t(temp)
    Sm <- t(sal)
    Depth = t(Depth)
  }else if (length(dim(temp))==3){
    n  <- dim(temp)[1]
    m  <- dim(temp)[2]
    d  <- dim(temp)[3]
    Tm <- t(matrix(temp,n*m,d))
    Sm <- t(matrix(sal,n*m,d))
  }else if(length(dim(temp))==4){
    temp <- aperm(temp, c(1,2,4,3))
    sal <- aperm(sal, c(1,2,4,3))
    n    <- dim(temp)[1]
    m    <- dim(temp)[2]
    t    <- dim(temp)[3]
    d    <- dim(temp)[4]
    Tm   <- t(matrix(temp,n*m*t,d))
    Sm   <- t(matrix(sal,n*m*t,d))
  }

  # first profile
  if(is.na(dim(Depth)[2])){
    depth <- Depth
    dm <- which.min(dmax-depth)
    k <- which(is.na(Tm[dm,])==FALSE)[1]
  }else{
    k <- 1
    depth <- depth[,k]
  }
  t <- Tm[,k]
  s <- Sm[,k]

  # Calcule les splines du profil 1
  myb     <- fda::create.bspline.basis(c(dmin,dmax),mybn)
  templ   <- approx(depth,t,method = 'linear',xout = profx)$y
  temp.fd <- fda::Data2fd(profx,templ,myb)
  sall    <- approx(depth,s,method = 'linear',xout = profx)$y
  sal.fd  <- fda::Data2fd(profx,sall,myb)
  ind <- k

  for (k in which(is.na(Tm[dm,])==FALSE)[2]:dim(Tm)[2]){  #In case of NaN in the matrix (land)
    t <- Tm[,k]
    s <- Sm[,k]
    if(is.na(t[dm])==FALSE){ #In case of NaN in the matrix (land)
      if(is.na(dim(Depth)[2])){
        depth = Depth
      }else{
        depth = Depth[,k]
      }
      templ         <- approx(depth,t,method = 'linear',xout = profx)$y
      tmpt.fd       <- fda::Data2fd(profx,templ,myb)
      sall          <- approx(depth,s,method = 'linear',xout = profx)$y
      tmps.fd       <- fda::Data2fd(profx,sall,myb)
      temp.fd$coefs <- cbind(temp.fd$coefs,tmpt.fd$coefs)
      sal.fd$coefs  <- cbind(sal.fd$coefs,tmps.fd$coefs)
      ind <- c(ind,k)
    }
  }

  Tfd <- fda::eval.fd(profx,temp.fd)
  Sfd <- fda::eval.fd(profx,sal.fd)

  #Creation of a mask of the map, when applicable
  if (length(dim(temp))==3){
    mask = matrix(NaN,n,m)
    mask[ind] = 1
  }else if(length(dim(temp))==4){
    mask = array(NaN,c(n,m,t))
    mask[ind] = 1
    mask = mask[,,1]
  }

  #return
  BS        <<- list(ind,profx,Tfd,Sfd)
  names(BS) <<- c('ind','profx','Teval','Seval')
  temp.fd   <<- temp.fd
  sal.fd    <<- sal.fd
  mask      <<- mask
  cat(paste(length(BS$ind)," profiles converted in Bsplines."))
}
