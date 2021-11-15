## Nov 1 


library(readxl)
Airbnb_Texas_Rentals <- read_excel("Airbnb_Texas_Rentals.xlsx", 
                                   col_types = c("numeric", "numeric", "numeric", 
                                                 "text", "text", "numeric", "numeric", 
                                                 "text", "numeric", "numeric", "text", 
                                                 "text"))

set.seed(100) 
airbnb=Airbnb_Texas_Rentals 
airbnb=airbnb[!is.na(airbnb$longitude) & !is.na(airbnb$latitude),] 
data=airbnb[airbnb$Year==2017,] 
data=data[!is.na(data$bedrooms_count),] 

library("fields") 
z=rnorm(2803,0,0.01) 
data$longitude=data$longitude+z 
# select a sample of 800
index=sample(1:2803, 800) 
#  select a sample the rest of the 2803, but not the 800 from index. 
# to have 2003 sample
index.v=setdiff(1:2803, index) 

data0=data[sample(1:2803, 800),] 
data.v=data[index.v,] # whats left over 
data.v=data.v[!is.na(data.v$bedrooms_count),] 


library("geoR") 

#fit=likfit(coords=cbind(data$longitude, data$latitude), data=log(data$average_rate_per_night), 
#           trend=trend.spatial(~data$bedrooms_count),ini.cov.pars=c(1,1), lik.method="ML") 


data1=data0[data0$longitude< -99,] 
data2=data0[data0$longitude >=-99 & data0$longitude< -97,] 
data3=data0[data0$longitude>= -97,] 


par(mfrow=c(1,1),mai=c(0.5,0.5,0.5,0.5)) 
quilt.plot(data0$longitude, data0$latitude, log(data0$average_rate_per_night)) 
US(add=T) 
abline(v=-99,col="gray") # ablines at -99 and -97 longatude
abline(v=-97, col="gray") 

text(-102, 28, label=nrow(data1), col="gray") 
text(-98, 28, label="418", col="gray") 
text(-95, 28, label=nrow(data3), col="gray") 

vario1=variog(coords=cbind(data1$longitude, data1$latitude), data=log(data1$average_rate_per_night), trend=trend.spatial(~data1$bedrooms_count),max.dist=3) 
vario2=variog(coords=cbind(data2$longitude, data2$latitude), data=log(data2$average_rate_per_night), trend=trend.spatial(~data2$bedrooms_count),max.dist=3) 
vario3=variog(coords=cbind(data3$longitude, data3$latitude), data=log(data3$average_rate_per_night), trend=trend.spatial(~data3$bedrooms_count),max.dist=3) 
par(mfrow=c(1,3),mai=c(0.5,0.5,0.5,0.5)) 
plot(vario1,ylim=c(0,1.1),main = "Data1") 
plot(vario2,ylim=c(0,1.1), main= "Data2") 
plot(vario3,ylim=c(0,1.1), main= "Data3") 

fit1=likfit(coords=cbind(data1$longitude, data1$latitude), data=log(data1$average_rate_per_night), 
            trend=trend.spatial(~data1$bedrooms_count),ini.cov.pars=c(1,1), lik.method="ML") 
fit2=likfit(coords=cbind(data2$longitude, data2$latitude), data=log(data2$average_rate_per_night), 
            trend=trend.spatial(~data2$bedrooms_count),ini.cov.pars=c(1,1), lik.method="ML") 
fit3=likfit(coords=cbind(data3$longitude, data3$latitude), data=log(data3$average_rate_per_night), 
            trend=trend.spatial(~data3$bedrooms_count),ini.cov.pars=c(1,1), lik.method="ML") 

print(list(fit1, fit2,fit3))

krig=function(data, data.v, par){ 
  n=dim(data)[1] 
  N=dim(data.v)[1] 
  alpha=par[1] 
  beta=par[2] 
  delta=par[3] 
  M=cbind(rep(1,n), data$bedrooms_count) 
  m=t(cbind(rep(1,N), data.v$bedrooms_count)) 
  D=rdist(cbind(data$longitude, data$latitude)) 
  d=rdist(cbind(data$longitude, data$latitude), cbind(data.v$longitude, data.v$latitude)) 
  S=alpha*exp(-D/beta) 
  diag(S)=diag(S)+delta 
  
  k=alpha*exp(-d/beta) 
  
  lambda=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
  krig=t(lambda) %*% log(data$average_rate_per_night) 
  
  return(krig) 
} 


## kriging with isotropic fit 
krig1=krig(data0, data.v,c(0.4339, 0.0552, 0.1999)) 



## kriging over region 1 
data.v1=data.v[data.v$longitude< -99,] 
krig2.1=krig(data1,data.v1, c(0.1878, 0.0624, 0.0324)) 

## kriging over region 2 
data.v2=data.v[data.v$longitude>= -99 & data.v$longitude< -97,] 
krig2.2=krig(data2,data.v2, c(0.4582, 0.0623, 0.0614)) 

## kriging over region 3 
data.v3=data.v[data.v$longitude>= -97,] 
krig2.3=krig(data3,data.v3, c(0.4930, 0.0399, 0.3025)) 



par(mfrow=c(1,1),mai=c(0.5,0.5,0.5,0.5)) 
plot(data0$longitude, data0$latitude, pch=20) 
US(add=T) 
abline(v=c(-99,-97),col="gray") 

points(data.v1$longitude, data.v1$latitude, col=2) 
points(data.v2$longitude, data.v2$latitude, col=3) 
points(data.v3$longitude, data.v3$latitude, col=4) 
par(mfrow=c(1,2),mai=c(0.5,0.8,0.5,0.8)) 

quilt.plot(data.v$longitude, data.v$latitude, log(data.v$average_rate_per_night) - krig1) 
US(add=T) 
data.vv=rbind(data.v1, data.v2, data.v3) 
krig.vv=c(krig2.1, krig2.2, krig2.3) 
quilt.plot(data.vv$longitude, data.vv$latitude, log(data.vv$average_rate_per_night) - krig.vv) 
US(add=T) 
par(mfrow=c(1,3),mai=c(0.5,0.5,0.5,0.5)) 
a=log(data.v$average_rate_per_night)-krig1 
b=log(data.vv$average_rate_per_night)-krig.vv 

boxplot(cbind(a,b)) 
plot(log(data.v$average_rate_per_night), krig1,pch=20,xlab="truth",ylab="krigged") 
points(log(data.vv$average_rate_per_night), krig.vv,col=2,pch=3) 
hist(a, freq=FALSE,ylim=c(0,2.5),main="",nclass=30,xlab="error") 
hist(b, freq=FALSE, add=T, border="red",nclass=30,col=NULL) 
