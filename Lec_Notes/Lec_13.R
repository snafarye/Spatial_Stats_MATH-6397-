## Oct 25, 2021 

# We will use Texas Airbnb data we used before for the exercise (see under “Data”) 
# Sample code that we will use as a basis is under “R code” 
# We will first explore appropriate spatial regression model 
# Then we will consider a few covariance models 
# We will use least squares method and likelihood based method to fit mean and covariance parameters 
# We will then perform Kriging on locations that we set aside 
# In the end, we can perform model comparison in terms of model fit as well as Kriging 
# We can also create Kriged map with associated uncertainty map 


# ----------------- Data clean up / exploration --------------------------------
library(readxl)
Airbnb_Texas_Rentals <- read_excel("Airbnb_Texas_Rentals.xlsx", 
                                   col_types = c("numeric", "numeric", "numeric", 
                                                 "text", "text", "numeric", "numeric", 
                                                 "text", "numeric", "numeric", "text", 
                                                 "text"))
airbnb = Airbnb_Texas_Rentals 

lapply(airbnb, function(x) { length(which(is.na(x)))}) # list of missing values in each column

airbnb=airbnb[!is.na(airbnb$longitude) & !is.na(airbnb$latitude),] # remove the NA values in Long and Lat columns
#hist(airbnb$average_rate_per_night)
log_price=log(airbnb$average_rate_per_night) 
#hist(log_price)

data=airbnb[airbnb$Year==2017,]  # just want to look at the year 2017

lapply(data, function(x) { length(which(is.na(x)))}) # see we have missing bedroom counts 
data=data[!is.na(data$bedrooms_count),] 
names(data)
dim(data) # 2803   12
#bedroom count, month, ...covaretae variables 

#-----------------------------------------------------------------
library("fields") 
set.seed(111)
# sample 500 rows 
index1=sample(1:dim(data)[1], 500) 
# 100 locations to do the kriging for prediction 
index2=sample(setdiff(1:dim(data)[1],index1), 100) 

data1=data[index1,] 
data2=data[index2,] 

# GOAL, dependent var = price , so what will be the x-var for our covariate??


# CHECK: is data normally distributed in order to use gaussian model? --> no its not, very skewed right
hist(data1$average_rate_per_night)  # to do the spacial regression model is the price normal? 
#------  to fit this, take the log, will be on the log scale
#   make sure to transform the data back taking the expo of price to have original scale of data
hist(log(data1$average_rate_per_night))
# much better distribution  


#----exploring spacial varb, price per night 
# variables to use for covariates variables ro help describle dependent variable ??
plot(data1$Month,          log(data1$average_rate_per_night)) # no clear linear/non-linear trend 
plot(data1$bedrooms_count, log(data1$average_rate_per_night)) # much more significant , pos trend 
plot(data1$longitude,      log(data1$average_rate_per_night)) # not really see a sig. relationship b/t price
plot(data1$latitude,       log(data1$average_rate_per_night)) # and long./ lat 
#---conclusion , pick bedroom counts and my only covareate and build a spacial cov model 
      #        UNIVERSAL  kriging since the mean is not constant over space 
#             y = int + bedroom_count(x)


fit0 = lm(log(data1$average_rate_per_night)~data1$bedrooms_count)
summary(fit0)
# any spatial dep in res? If not --> no need spacial modeling. 
# If yes? then need to further incorpurate spatial dep in model
res0 = fit0$residuals  
# spatial map, see any spatial patterns? Hard to see from map
quilt.plot(data1$longitude,
           data1$latitude, 
           res0)
library("geoR")
vario0 = variog(coords = cbind(data1$longitude, data1$latitude), 
                data = res0)
plot(vario0) # have very short spatial dependence 
# euclidean distance using degrees




 
#----------  universal kriging  ----------------------------------------------

## variogram  (original data, not the res) 
vario1=variog(coords=cbind(data1$longitude, data1$latitude), 
              data=log(data1$average_rate_per_night)) 
plot(vario1) 

vario2=variog(coords=cbind(data1$longitude, data1$latitude),
              data=log(data1$average_rate_per_night), 
              trend=trend.spatial(~data1$bedrooms_count), 
              max.dist=6) 
plot(vario2) 



##---------- least squares ----------------------------------------------------
# geoR package 
# smaller minimized = better 
fit.expo=variofit(vario2, ini.cov.pars=c(1,1)) ## try many other models 
fit.expo
#  tausq  sigmasq     phi 
#  0.0000  0.6792     0.1160 
#  Practical Range with cor=0.05 for asymptotic range: 0.347542
#  variofit: minimised weighted sum of squares = 300.5722


fit.shp=variofit(vario2, ini.cov.pars=c(1,1),
                 cov.model = "sph") ## try many other models 
fit.shp
# expo is better since it have smaller min weighted sum of squares 
#  tausq   sigmasq     phi 
#  0.0000  0.6793      0.3410 
#  Practical Range with cor=0.05 for asymptotic range: 0.3409595
#  variofit: minimised weighted sum of squares = 298.8389



##---------  likelihood method ------------------------------------------------
# MLE vs REML 

#ML, expo
fit2=likfit(coords=cbind(data1$longitude, data1$latitude), 
            data=log(data1$average_rate_per_night), 
            trend=trend.spatial(~data1$bedrooms_count), # mean structure
            lik.method = "ML",
            ini.cov.pars=c(0.5,2))  
# error due to duplicates in the months of year 2017, multiple measurements 
# can not calc inverse , get a large nugget 
# now how to remove the duplicates 

##-------- duplicates! 
zz=rnorm(dim(data1)[1], 0,0.1) 
data1$longitude=data1$longitude+zz 
#run again
fit2=likfit(coords=cbind(data1$longitude, data1$latitude), 
            data=log(data1$average_rate_per_night), 
            trend=trend.spatial(~data1$bedrooms_count), # mean structure
            lik.method = "ML",
            ini.cov.pars=c(0.5,2)) 
fit2# now it works 


# MLe, spherical  
fit3=likfit(coords=cbind(data1$longitude,data1$latitude), 
            data=log(data1$average_rate_per_night), 
            trend=trend.spatial(~data1$bedrooms_count),
            cov.model="spherical", 
            ini.cov.pars=c(0.5,2) ) 
fit3

# expo is better with a larger max log-likelihood 




#------------------Example reml--------------
fit2.reml=likfit(coords=cbind(data1$longitude, data1$latitude), 
            data=log(data1$average_rate_per_night), 
            trend=trend.spatial(~data1$bedrooms_count), # mean structure
            lik.method = "REML",
            ini.cov.pars=c(0.5,2)) 
fit2.reml# no need to do that duplicates part like above for reml

# MLe, spherical  
fit3.reml=likfit(coords=cbind(data1$longitude,data1$latitude), 
            data=log(data1$average_rate_per_night), 
            trend=trend.spatial(~data1$bedrooms_count),
            lik.method = "REML",
            cov.model="spherical", 
            ini.cov.pars=c(0.5,2) ) 
fit3.reml


# reml is also better since is it unbiased for mle estimates 
# again see that expo is better with a larger max log-likelihood 




##---------- Kriging  --------------------------------------------------------- 
#   remeber      y = int + bedroom_count(x)
# modify code below  
M=cbind(rep(1,500),  data1$bedrooms_count) 
m=t(cbind(rep(1,100), data2$bedrooms_count)) 

#reml results with exponential cov fun (prof results)
#alpha = 0.1744
#beta = 0.6234
#delta = 0.4333
# -------------------------(my results)
alpha = 0.2680
beta = 0.2219
delta = 0.3786

D = rdist(cbind(data1$longitude,data1$latitude))   # distance matrix 
S = alpha * exp(-D/beta)                           # cov matrix 
diag(S) = diag(S)+delta                            # add nugget 

d = rdist(cbind(data1$longitude, data1$latitude),  # to calc little k
          cbind(data2$longitude, data2$latitude))
k = alpha*exp(-d/beta)
lambda = (solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 


krig.expo.reml = t(lambda) %*% log(data1$average_rate_per_night) 
# predict values in log
plot(log(data2$average_rate_per_night), krig.expo.reml) #x --> true log price, y --> kriging results, 45 degree  log scale

mse = mean((krig.expo.reml-log(data2$average_rate_per_night))^2);mse
mae = mean(abs(krig.expo.reml-log(data2$average_rate_per_night)));mae



## Krigged map with fields package 
# 
fit=Krig(cbind(data1$longitude, data1$latitude), log(data1$average_rate_per_night),m=1, Z=data1$bedrooms_count) 
surface(fit) 
out.p1=predict(fit, cbind(data2$longitude, data2$latitude),Z=data2$bedrooms_count) 
out.p2=predictSE(fit, cbind(data2$longitude, data2$latitude),Z=data2$bedrooms_count) 
quilt.plot(cbind(data2$longitude, data2$latitude), out.p1) 
quilt.plot(cbind(data2$longitude, data2$latitude), out.p2) 

data.all=data[setdiff(1:dim(data)[1], index1),] 
out.p3=predict(fit, cbind(data.all$longitude, data.all$latitude),Z=data.all$bedrooms_count) 
out.p4=predictSE(fit, cbind(data.all$longitude, data.all$latitude),Z=data.all$bedrooms_count) 
quilt.plot(cbind(data.all$longitude, data.all$latitude), out.p3-log(data.all$average_rate_per_night)) 







## =============================================================================

Library(sp) 

##### REML code example ##### 
# gudicdist, distance matrix

D=rdist.earth(cbind(data1$longitude, data1$latitude),miles=F) 
m=2 # number of mean parameters 

reml.loglik=function(par){ 
  print(par) 
  alpha=exp(par[1])
  beta=exp(par[2]) 
  nu= exp(par[3])/(1+exp(par[3]))*5 
  delta=exp(par[4])   
  
  M=cbind(rep(1, dim(D)[1]), data1$bedrooms_count)    ## design matrix 
  Z=matrix(log(data1$average_rate_per_night),ncol=1)  # spatial data vector 
  S=Matern(D/beta, smoothness=nu)*alpha               # matren cov matrix 
  diag(S)=diag(S)+delta                               # nugget 

  C=log(2*pi)*(-(dim(D)[1]-m)/2)                      # calc log reml likelihood
  Y=( diag(1, dim(D)[1],dim(D)[1])-  M %*% solve(t(M) %*% M) %*% t(M) ) %*% Z 
  
  t1=determinant(S, logarithm=T)$modulus 
  t2=determinant( t(M) %*% solve(S) %*% M, logarithm=T)$modulus 
  
  t3=t(Y) %*% (solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% Y 
  temp= -C - (t1+t2+t3)/2 
  print(temp) 
  return(-temp) 
} 

ini=c(0.1,5,-1,-3) 
fit=optim(ini, reml.loglik,control=list(maxit=100000000)) 

#> fit 

#$par 
#[1] -1.5215103  3.9685613 -3.1878741 -0.8559631 

#$value 
#[1] -338.851 


#$counts 
#function gradient  
#327       NA  

#$convergence 
#[1] 0 

#$message 
#NULL 


