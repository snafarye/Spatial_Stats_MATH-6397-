
#################### START OF SPHERICAL OLS ON LOCAL STATIONARY  #########################
##-----------  local statinarity ----------------------------

# quick renaming 
data0  = sample_data
data.v = sample_kriging_data

range(sample_data$longitude)

data1=data0[data0$longitude > -122.204,] 
data2=data0[data0$longitude <= -122.204 & data0$longitude >  -122.397,]
data3=data0[data0$longitude <=  -122.397 & data0$longitude> -122.59,]

par(mfrow=c(1,1)) 
quilt.plot(data0$longitude, data0$latitude, res0) 
US(add=T) 
abline(v= c( -122.204, -122.397),col="gray") # ablines at -122.204, -122.397 longatude


text( c(-122.5, -122.3, -122.15), 
      c(37.6, 37.6, 37.6), 
      label = c(nrow(data3), nrow(data2), nrow(data1) ), 
      col = "black")

text( c(-122.5, -122.3, -122.15), 
      c(37.65, 37.65, 37.65), 
      label = c("data3", "data2", "data1" ), 
      col = "black")


print(
  cat(" Number of cases within each subgroup:", '\n',
      "data1 = ", nrow(data1) ,'\n',
      "data2 = ", nrow(data2) ,'\n',
      "data3 = ", nrow(data3) ,'\n'
  )
)

###-------- variograms of the subgroups-----------------------------
res_krig_data1 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data1))$residuals

vario1=variog(coords=cbind(data1$longitude, data1$latitude), 
              data=res_krig_data1, 
              trend=trend.spatial(~data1$total_bedrooms  + 
                                    data1$median_income +
                                    data1$housing_median_age +
                                    data1$longitude+ 
                                    data1$latitude
              ),
              max.dist = 0.6) 
plot(vario1)

#SPHERICAL 
plot(vario1, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.08,0.2), 
                 nugget = 0.005, 
                 max.dist = 0.6,
                 lwd = 3)

res_krig_data2 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data2))$residuals

vario2=variog(coords=cbind(data2$longitude, data2$latitude), 
              data=res_krig_data2, 
              trend=trend.spatial(~data2$total_bedrooms  + 
                                    data2$median_income +
                                    data2$housing_median_age +
                                    data2$longitude+ 
                                    data2$latitude),
              max.dist = 0.6)
plot(vario2)
plot(vario2, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.1,0.1), 
                 nugget = 0.005, 
                 max.dist = 1,
                 lwd = 3)


res_krig_data3 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data3))$residuals

vario3=variog(coords=cbind(data3$longitude, data3$latitude), 
              data=res_krig_data3, 
              max.dist = 0.39,
              trend=trend.spatial(~data3$total_bedrooms  + 
                                    data3$median_income +
                                    data3$housing_median_age +
                                    data3$longitude+ 
                                    data3$latitude
              ))
plot(vario3)
plot(vario3, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.08,0.10), 
                 nugget = 0.005, 
                 max.dist = 0.39,
                 lwd = 3)




#---plot-----------------------------------------------------

par(mfrow=c(2,2)) 
plot(vario1, main = "Data1")
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.08,0.2), 
                 nugget = 0.0005, 
                 max.dist = 0.6,
                 lwd = 3)

plot(vario2, main= "Data2") 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.10,0.1), 
                 nugget = 0.0005, 
                 max.dist = 0.6,
                 lwd = 3)
plot(vario3, main= "Data3")
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.08,0.10), 
                 nugget = 0.0005, 
                 max.dist = 0.39,
                 lwd = 3)

par(mfrow=c(1,1))

###------ SPHERICAL OLS  variofit() of the subgroups-------------------------------------

fit1=variofit(vario1, ini.cov.pars=c(0.08,0.2),
              weights = 'equal',
              cov.model = "spherical")
fit1
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0132  0.0272  0.0458 
#Practical Range with cor=0.05 for asymptotic range: 0.0458001

#variofit: minimised sum of squares = 0.0019

fit2=variofit(vario2, ini.cov.pars=c(0.1,0.1),
              weights = 'equal',
              cov.model = "spherical")
fit2

#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0682  0.0430  0.0621 
#Practical Range with cor=0.05 for asymptotic range: 0.0620967

#variofit: minimised sum of squares = 0.0223



fit3=variofit(vario3, ini.cov.pars=c(0.08,0.1),
              weights = 'equal',
              cov.model = "spherical")
fit3  

#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0340  0.0619  0.0324 
#Practical Range with cor=0.05 for asymptotic range: 0.03235643

#variofit: minimised sum of squares = 0.0012


print(list(fit1, fit2, fit3))
> print(list(fit1, fit2, fit3))
[[1]]
variofit: model parameters estimated by OLS (ordinary least squares):
  covariance model is: spherical
parameter estimates:
  tausq sigmasq     phi 
0.0368  0.0468  0.2288 
Practical Range with cor=0.05 for asymptotic range: 0.2287894

variofit: minimised sum of squares = 9e-04

[[2]]
variofit: model parameters estimated by OLS (ordinary least squares):
  covariance model is: spherical
parameter estimates:
  tausq sigmasq     phi 
0.0352  0.0734  0.0443 
Practical Range with cor=0.05 for asymptotic range: 0.04431049

variofit: minimised sum of squares = 0.0021

[[3]]
variofit: model parameters estimated by OLS (ordinary least squares):
  covariance model is: spherical
parameter estimates:
  tausq sigmasq     phi 
0.0557  0.0253  0.0335 
Practical Range with cor=0.05 for asymptotic range: 0.03348256

variofit: minimised sum of squares = 0.0021

#------- SPHERICAL OLS KRIGING FOR LOCAL STATIONARY

# quick reminder 
data0  = sample_data
data.v = sample_kriging_data

#-- function

krig=function(data, data.v, par){ 
  n=dim(data)[1] 
  N=dim(data.v)[1] 
  alpha=par[1] 
  beta=par[2] 
  delta=par[3] 
  M=cbind(rep(1,n), data$total_bedrooms,data$median_income,data$housing_median_age,data$longitude,data$latitude) 
  m=t(cbind(rep(1,N), data.v$total_bedrooms,data.v$median_income,data.v$housing_median_age,data.v$longitude,data.v$latitude)) 
  D=rdist(cbind(data$longitude, data$latitude)) 
  d=rdist(cbind(data$longitude, data$latitude), cbind(data.v$longitude, data.v$latitude)) 
  S=alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3)))) 
  diag(S)=diag(S)+delta 
  
  
  
  k=alpha*(1-((3*d)/(2*beta))+((d**3)/(2*(beta**3))))
  T = 0
  lambda = (solve(S, tol = T) - solve(S, tol =T) %*% M %*% solve(t(M) %*% solve(S, tol =T) %*% M,tol = T) %*% t(M) %*% solve(S,tol = T) ) %*% k + solve(S,tol = T) %*% M %*% solve(t(M) %*% solve(S,tol =T) %*% M,tol =T) %*% m 

  krig=t(lambda) %*% lm(log(median_house_value) ~ total_bedrooms+ 
                          median_income+housing_median_age + 
                          longitude+latitude,  data = data)$residuals
  
  return(krig) 
} 


# ----   results 


## kriging with isotropic fit 
krig0 = krig(data0, data.v, c( 0.0587, 0.0520, 0.0480)) 

#--------------------------------------------------------change the range 

print(list(fit1, fit2, fit3))

## kriging over region 1 
data.v1 = data.v[data.v$longitude > -122.204,] 
krig2.1 = krig(data1, data.v1, c(0.0468  ,  0.2288  , 0.0368  )) 

## kriging over region 2 
data.v2 = data.v[data.v$longitude <= -122.204 & data.v$longitude >  -122.397,] 
krig2.2 = krig(data2, data.v2, c(0.0734  ,  0.0443 , 0.0352    )) 

## kriging over region 3 
data.v3 = data.v[data.v$longitude <=  -122.397 & data.v$longitude > -122.59,] 
krig2.3 = krig(data3, data.v3, c(0.0253  ,  0.0335 ,0.0557   ))


error = function(data, krig_name){
  res_krig <- (lm(log(median_house_value) ~ total_bedrooms+ 
                    median_income+housing_median_age + 
                    longitude+latitude,  data = data))$residuals
  res_krig = matrix(res_krig, ncol = 1)
  plot(res_krig, krig_name, xlab = "true", ylab = "kriged value")
  lines(res_krig,res_krig) #y=x line 
  
  mse=mean((krig_name-res_krig)^2)
  mae=mean(abs(krig_name-res_krig))
  
  list_error <- list("MSE" = mse, "MAE" = mae)
  return(list_error)
  
  
}

error(data.v1, krig2.1)
error(data.v2, krig2.2)
error(data.v3, krig2.3)

> error(data.v1, krig2.1)
$MSE
[1] 0.03122966

$MAE
[1] 0.1316569

> error(data.v2, krig2.2)
$MSE
[1] 0.03941714

$MAE
[1] 0.1496382

> error(data.v3, krig2.3)
$MSE
[1] 0.02333017

$MAE
[1] 0.1189648



# plotting the kriging 
par(mfrow=c(1,1)) 
plot(data0$longitude, data0$latitude, pch=20, main = "SPHERICAL OLS") 
US(add=T) 
abline(v=c(-122.397, -122.204),col="gray") 

points(data.v1$longitude, data.v1$latitude, col=2, lwd = 3) 
points(data.v2$longitude, data.v2$latitude, col=3, lwd = 3) 
points(data.v3$longitude, data.v3$latitude, col=4, lwd = 3) 



# --------checking errors v true values----------------------------- 
par(mfrow=c(1,2))

quilt.plot(data.v$longitude, data.v$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                   median_income+housing_median_age + 
                                                   longitude+latitude,  data = data.v)$residuals - krig0,
           main= "overall") 
US(add=T) 
data.vv=rbind(data.v1, data.v2, data.v3) 
krig.vv=c(krig2.1, krig2.2, krig2.3) 
quilt.plot(data.vv$longitude, data.vv$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                     median_income+housing_median_age + 
                                                     longitude+latitude,  data = data.vv)$residuals - krig.vv,
           main = "local stationarity") 
US(add=T)


par(mfrow=c(1,3),mai=c(0.5,0.5,0.5,0.5)) 
a=lm(log(median_house_value) ~ total_bedrooms+ 
       median_income+housing_median_age + 
       longitude+latitude,  data = data.v)$residuals - krig0 
b=lm(log(median_house_value) ~ total_bedrooms+ 
       median_income+housing_median_age + 
       longitude+latitude,  data = data.vv)$residuals-krig.vv 


boxplot(cbind(a,b)) 
plot(lm(log(median_house_value) ~ total_bedrooms+ 
          median_income+housing_median_age + 
          longitude+latitude,  data = data.v)$residuals , krig0,pch=20,xlab="truth",ylab="krigged") 
points(lm(log(median_house_value) ~ total_bedrooms+ 
            median_income+housing_median_age + 
            longitude+latitude,  data = data.vv)$residuals, krig.vv,col=2,pch=3, lwd = 2) 
hist(a, freq=FALSE,main="",nclass=30,xlab="error") 
hist(b, freq=FALSE, add=T, border="red",nclass=30,col=NULL, lwd = 2) 



#################### END OF SPHERICAL OLS ON LOCAL STATIONARY  #########################


#################### START OF SPHERICAL WLS ON LOCAL STATIONARY  #########################


###------ SPHERICAL WLS  variofit() of the subgroups-------------------------------------

fit1a=variofit(vario1, ini.cov.pars=c(0.08,0.2),
               weights = 'npairs',
               cov.model = "spherical")
fit1a
#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0133  0.0321  0.0558 
#Practical Range with cor=0.05 for asymptotic range: 0.05582116

#variofit: minimised weighted sum of squares = 1.4928

fit2a=variofit(vario2, ini.cov.pars=c(0.1,0.1),
               weights = 'npairs',
               cov.model = "spherical")
fit2a

#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0678  0.0591  0.0873 
#Practical Range with cor=0.05 for asymptotic range: 0.08727192

#variofit: minimised weighted sum of squares = 70.0209



fit3a=variofit(vario3, ini.cov.pars=c(0.08,0.1),
               weights = 'npairs',
               cov.model = "spherical")
fit3a  

#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0342  0.0633  0.0334 
#Practical Range with cor=0.05 for asymptotic range: 0.03337693

#variofit: minimised weighted sum of squares = 9.7941


#------- SPHERICAL WLS KRIGING FOR LOCAL STATIONARY

#-- function


# ----   results 

print(list(fit1a,fit2a,fit3a))
> print(list(fit1a,fit2a,fit3a))
[[1]]
variofit: model parameters estimated by WLS (weighted least squares):
  covariance model is: spherical
parameter estimates:
  tausq sigmasq     phi 
0.0537  0.0323  0.4306 
Practical Range with cor=0.05 for asymptotic range: 0.4305687

variofit: minimised weighted sum of squares = 4.3566

[[2]]
variofit: model parameters estimated by WLS (weighted least squares):
  covariance model is: spherical
parameter estimates:
  tausq sigmasq     phi 
0.0354  0.0748  0.0455 
Practical Range with cor=0.05 for asymptotic range: 0.04552278

variofit: minimised weighted sum of squares = 15.7368

[[3]]
variofit: model parameters estimated by WLS (weighted least squares):
  covariance model is: spherical
parameter estimates:
  tausq sigmasq     phi 
0.0646  0.0186  0.0769 
Practical Range with cor=0.05 for asymptotic range: 0.07686287

variofit: minimised weighted sum of squares = 7.5998

## kriging with isotropic fit 

krig0a = krig(data0, data.v, c( 0.0593, 0.0663, 0.0480)) 

#--------------------------------------------------------change the range 

print(list(fit1a, fit2a, fit3a))

## kriging over region 1 
data.v1 = data.v[data.v$longitude > -122.204,] 
krig2.1 = krig(data1,data.v1, c(0.0323    ,  0.4306  , 0.0537    )) 

## kriging over region 2 
data.v2 = data.v[data.v$longitude <= -122.204 & data.v$longitude>  -122.397,] 
krig2.2 = krig(data2,data.v2, c(0.0748    ,  0.0455  , 0.0354      )) 

## kriging over region 3 
data.v3 = data.v[data.v$longitude <=  -122.397 & data.v$longitude > -122.59,] 
krig2.3 = krig(data3,data.v3, c(0.0253  ,  0.0335 ,0.0557   ))



error = function(data, krig_name){
  res_krig <- (lm(log(median_house_value) ~ total_bedrooms+ 
                    median_income+housing_median_age + 
                    longitude+latitude,  data = data))$residuals
  res_krig = matrix(res_krig, ncol = 1)
  plot(res_krig, krig_name, xlab = "true", ylab = "kriged value")
  lines(res_krig,res_krig) #y=x line 
  
  mse=mean((krig_name-res_krig)^2)
  mae=mean(abs(krig_name-res_krig))
  
  list_error <- list("MSE" = mse, "MAE" = mae)
  return(list_error)
  
  
}

error(data.v1, krig2.1)
error(data.v2, krig2.2)
error(data.v3, krig2.3)

> error(data.v1, krig2.1)
$MSE
[1] 0.03422007

$MAE
[1] 0.1404673

> error(data.v2, krig2.2)
$MSE
[1] 0.03923733

$MAE
[1] 0.1494093

> error(data.v3, krig2.3)
$MSE
[1] 0.02333017

$MAE
[1] 0.1189648




# plotting the kriging 
par(mfrow=c(1,1)) 
plot(data0$longitude, data0$latitude, pch=20, main = "SPHERICAL WLS") 
US(add=T) 
abline(v=c(-122.204,-122.397, -122.59),col="gray") 

points(data.v1$longitude, data.v1$latitude, col=2, lwd = 3) 
points(data.v2$longitude, data.v2$latitude, col=3, lwd = 3) 
points(data.v3$longitude, data.v3$latitude, col=4, lwd = 3) 

# --------checking errors v true values----------------------------- 
par(mfrow=c(1,2))

quilt.plot(data.v$longitude, data.v$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                   median_income+housing_median_age + 
                                                   longitude+latitude,  data = data.v)$residuals - krig0,
           main= "overall") 
US(add=T) 
data.vv=rbind(data.v1, data.v2, data.v3) 
krig.vva=c(krig2.1a, krig2.2a, krig2.3a) 
quilt.plot(data.vv$longitude, data.vv$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                     median_income+housing_median_age + 
                                                     longitude+latitude,  data = data.vv)$residuals - krig.vv,
           main = "local stationarity") 
US(add=T)


par(mfrow=c(1,3),mai=c(0.5,0.5,0.5,0.5)) 
a=lm(log(median_house_value) ~ total_bedrooms+ 
       median_income+housing_median_age + 
       longitude+latitude,  data = data.v)$residuals - krig0a 
b=lm(log(median_house_value) ~ total_bedrooms+ 
       median_income+housing_median_age + 
       longitude+latitude,  data = data.vv)$residuals-krig.vva 


boxplot(cbind(a,b)) 
plot(lm(log(median_house_value) ~ total_bedrooms+ 
          median_income+housing_median_age + 
          longitude+latitude,  data = data.v)$residuals , krig0a,pch=20,xlab="truth",ylab="krigged") 
points(lm(log(median_house_value) ~ total_bedrooms+ 
            median_income+housing_median_age + 
            longitude+latitude,  data = data.vv)$residuals, krig.vva,col=2,pch=3, lwd = 2) 
hist(a, freq=FALSE,main="",nclass=30,xlab="error") 
hist(b, freq=FALSE, add=T, border="red",nclass=30,col=NULL, lwd = 2) 



#################### END OF SPHERICAL OLS ON LOCAL STATIONARY  #########################



