library(devtools)
devtools::install_github("Epauthenet/fda.oce")
#help(package = fda.oce)
library(fda.oce)
###########################

install.packages("roxygen2")
install.packages("devtools")
library(roxygen2)
library(devtools)
has_devel()

devtools::document("~/Documents/R/fda.oce")

#
#DONE
? bspl
? fpca
? reco
? pc_plot
? eigenf_plot
? reco
#
.libPaths()
#

#####Test run
setwd("~/Documents/R/fda.oce")
library(ramify)
load("GLORYS_2015-12_SO_sub10.RData")

Temp = t(resize(temp,144*17,47))
Sal = t(resize(sal,144*17,47))

Xi = array(c(Temp,Sal),c(47,144*17,2))
Xi = Xi[,colSums(is.na(Xi[,,1]))<47,]   #Remove lands
Xi = Xi[,is.na(Xi[47,,1])==F,]   #Remove short profiles
Pi = as.numeric(depth)
save(Pi,Xi,file = "GLORYS_SO_2015-12.RData")
library(R.matlab)
writeMat("GLORYS_SO_2015-12.mat",Pi = Pi,Xi = Xi)
#
#####################################
load("GLORYS_SO_2015-12.RData")
fda.oce::bspl(Pi,Xi)

#
fpca(fdobj)
proj(fdobj,pca)

png("~/Documents/R/fda.oce/figures/pc_plot.png", width=7, height=7, units="in", res=600)
pc_plot(pca,pc,c(1,2))
dev.off()
#
for(te in 1:2){
  png(paste("~/Documents/R/fda.oce/figures/eigenf",te,".png",sep = ""), width=7, height=7, units="in", res=600)
  eigenf_plot(pca,te)
  dev.off()
}

te = 5
reco(pca,pc,te)

data = eval.fd(Pi,fdobj)
data_reco = eval.fd(Pi,fdobj_reco)

i = 600  #index of a profile
png(paste("~/Documents/R/fda.oce/figures/reco_prof",i,".png",sep = ""), width=7, height=7, units="in", res=600)
par(mfrow = c(1,2))
for(k in 1:pca$ndim){
  plot(Xi[,i,k],Pi,las = 1,cex = .2,col = 1
    ,xlim = range(Xi[,i,k],data_reco[,i,k])
    ,ylim = c(1000,0)
    ,xlab = pca$fdnames[[2+k]]
    ,ylab = pca$fdnames[[1]])
  points(data[,i,k],Pi,typ = 'l',col = 2)
  points(data_reco[,i,k],Pi,las = 1,ylim = c(1000,0),typ = 'l',col = 3)
}
legend("bottomleft",col = c(1,2,3),lty = c(NA,1,1),pch = c(20,NA,NA),cex = .7
  ,legend = c("raw data","B-spline fit",paste("reconstruction with ",te," modes",sep = "")))
dev.off()
#










##########The test run will be done on a subsample of one monthly mean field from Global Ocean Physics Reanalysis
#(mercatorglorys2v4_gl4_mean_201512.nc)

library(ncdf4)
setwd("~/Documents/Database/GLORYS_2015")
#Load the netcdf
ncid  <- nc_open("mercatorglorys2v4_gl4_mean_201512.nc")
Sal   <- ncvar_get( ncid,'salinity')
Temp  <- ncvar_get( ncid,'temperature')-273.15
Depth <- ncvar_get( ncid,'depth')
Lat   <- ncvar_get( ncid,'latitude')
Lon   <- ncvar_get( ncid,'longitude')

#Take a subsample
lo <- seq(1,length(Lon),10)
la <- seq(1,161,10)
de <- 1:47
temp  <- Temp[lo,la,de]
sal   <- Sal[lo,la,de]
lon   <- Lon[lo]
lat   <- Lat[la]
depth <- Depth[de]
rm(temp,sal)
nc_close(ncid)

save(temp,sal,lon,lat,depth,ncid,mask,file = "GLORYS_2015-12_SO_sub10.RData")
#---------------------------------------------------------------------------------------Bsplines
setwd("~/Documents/R/fda.oce/")
load("~/Documents/Database/GLORYS_2015/GLORYS_2015-12_SO_sub10.RData")

fda.oce::bspl(temp,sal,depth,range = c(5,1000),mybn = 20)

load("~/Documents/Database/GLORYS_2015/GLORYS_2015-12_SO_sub10_BS.RData")

fda.oce::fpca(temp.fd,sal.fd,plot = T)

reco(pca)

load("i06s_1000m40BS.RData")
proj(ti06.fd,si06.fd,pca)


png("~/Documents/R/fda.oce/GLO_PCmap.png", width=12, height=5, units="in", res=600)
par(mfrow = c(2,2),mar = c(1,3,3,3))
for (te in 1:4){
  map_pc(lon,lat,pca,mask,te,sign = c(1,1,1,1))
}
dev.off()

png("~/Documents/R/fda.oce/GLO_kde.png", width=8, height=4, units="in", res=600)
par(mfrow = c(1,2),mar = c(3,3,3,2))
plot(pca$pc[,1:2],xlab=paste("PC1 (",pca$pval[1]," %)"),las = 1,main = "Southern Ocean (GLORYS) 12-2015"
  ,ylab=paste("PC2 (",pca$pval[2]," %)"),col = 1,pch = 20,cex = .4)
grid()
abline(v=0,h=0,lty=3)
#
par(mar = c(3,1,3,5))
kde_pc(pca, te = c(1, 2), h = 0.2, n = 50,title = "KDE and projection")
points(Npc[,1:2],col = 2,pch = "+")
legend("topright",legend = "section WOCE i06s",pch = "+",col = 2)
dev.off()
#

for (te in 1:4){
  png(paste("~/Documents/R/fda.oce/GLO_eigen",te,".png",sep = ""), width=5, height=5, units="in", res=600)
  eigenf(pca,te)
  dev.off()
}
#

#############################


save(lon,lat,depth,mask,ncid,BS,temp.fd,sal.fd,file = "GLORYS_2015-12_SO_sub10_BS.RData")
save(lon,lat,depth,mask,pca,file = "GLORYS_2015-12_SO_sub10_pca.RData")


















#######################Preparation of i06s
library(ncdf4)
library(maps)
library(oce)
library(fields)
library(fda)
#

#
lat = NULL
lon = NULL
time = list()
profx = seq(5,1000)  #only down to 1000m
Temp = matrix(NaN,996,213)
Psal = matrix(NaN,996,213)
Depth = matrix(NaN,996,213)
lf = list.files("~/Documents/Database/WOCE/i06s_33RR20080204_nc_ctd",pattern = ".nc") #(Pacific)
t = 1
for (t in 1:length(lf)){
  ncid = nc_open(paste("~/Documents/Database/WOCE/i06s_33RR20080204_nc_ctd/",lf[t],sep = ""))
  depth = ncvar_get(ncid,'pressure')
  Temp[,t] = approx(depth,ncvar_get(ncid,'temperature'),profx)$y
  Psal[,t] = approx(depth,ncvar_get(ncid,'salinity'),profx)$y
  lat[t] = ncvar_get(ncid,'latitude')
  lon[t] = ncvar_get(ncid,'longitude')
  time[[t]] = as.Date(round(ncvar_get(ncid,'time')/1440,0),origin = '1980-01-01') #units: minutes since 1980-01-01 00:00:00
  nc_close(ncid)
}

#
i = 25:78
plot(lon[i],lat[i],xlim = c(0,50))
map(add = T)

#####################################Bsplines
##############add data points at depth for one profile
Temp[947:996,i[70]] = Temp[946,i[70]]
Psal[947:996,i[70]] = Psal[946,i[70]]
dmin = 5
dmax = 1000
mybn = 20
myb  = create.bspline.basis(c(dmin,dmax),mybn)
profx = seq(dmin,dmax)
k = i[1]
temp.fd = Data2fd(profx,Temp[,k],myb)
sal.fd = Data2fd(profx,Psal[,k],myb)
for (k in i[-1]){
  tmpt.fd = Data2fd(profx,Temp[,k],myb)
  tmps.fd  = Data2fd(profx,Psal[,k],myb)
  temp.fd$coefs = cbind(temp.fd$coefs,tmpt.fd$coefs)
  sal.fd$coefs  = cbind(sal.fd$coefs,tmps.fd$coefs)
  print(k)
}
Tfd = eval.fd(profx,temp.fd)
Sfd = eval.fd(profx,sal.fd)
#save(ti06.fd,si06.fd,file = "~/Documents/R/fda.oce/i06s_1000m40BS.RData")
