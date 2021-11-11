library(readxl)
load("C:/Users/Samue/Desktop/books for school/FALL 2021/MATH 6397/data1.RData")
dim(data1)
library(fields)
#library(gstat) 
#library(sp) 
library(geoR)

#---------------------- Data Cleaning ------------------------------------------

general_data = data1 #Create a dataframe with the loaded dataset
nrow(general_data) #Check the number of rows of the dataset
sum(is.na(general_data)) # have about 20 rows with at least one missing value
colSums(is.na(general_data))  # check to see what columns have the missing values 

new_data <- na.omit(general_data) #create a new dataset without missing values
nrow(new_data) #check number of rows of new dataset without missing values
sum(is.na(new_data)) #check to see if all missing values have been removed

set.seed(222)
# 1500 sample data
sample_data = new_data[sample(1:dim(new_data)[1],1500), ] #create a sample dataset of 1500 rows
nrow(sample_data)#check to see if sample dataset has 1500 rows
# 100 kriging data 
sample_kriging_data = new_data[sample(setdiff(1:dim(new_data)[1], sample_data),100),] # 100 location for kriging
nrow(sample_kriging_data)

loc = cbind(sample_data$longitude, sample_data$latitude)
y      <- log(sample_data$median_house_value) #value at loc
x1     <- sample_data$total_bedrooms 
x2     <- sample_data$median_income  
x3     <- sample_data$housing_median_age 
x4     <- sample_data$longitude
x5     <- sample_data$latitude

# linear regression model from report 1
fit0 = lm(log(median_house_value) ~ total_bedrooms+ 
            median_income+housing_median_age + 
            longitude+latitude,  data = sample_data)
summary(fit0)
res0 = fit0$residuals

# spatial map
quilt.plot(sample_data$longitude, sample_data$latitude, y)
US(add=T) 


#----- linear regression model from report 1
fit0 = lm(y ~ x1+x2+x3+x4+x5,  data = sample_data)
summary(fit0)
res0 = fit0$residuals
summary(res0)

plot(res0)
quilt.plot(sample_data$longitude, sample_data$latitude, res0)
US(add=T) 


##########  START OF  SPHERICAL  FOR OVERALL DATA  #########################

#------------------- variogram----------------
vario0  = variog(coords=loc,
                 data=res0, 
                 trend=trend.spatial(~x1+x2+x3+x4+x5),
                 max.dist = 0.6) 
plot(vario0) 




# spherical 
plot(vario0, main = "SHPERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.1,0.05), 
                 nugget = 0.01, 
                 max.dist = 0.6,
                 lwd = 3)

#################### START OF SPHERICAL OLS ON OVERALL  #########################

fit.sph.OLS=variofit(vario0, ini.cov.pars=c(0.15,0.05),
                      weights = 'equal',
                      cov.model = "spherical")
fit.sph.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0480  0.0587  0.0520 
#Practical Range with cor=0.05 for asymptotic range: 0.05196297

#variofit: minimised sum of squares = 0.0016

set.seed(222)
zz=rnorm(dim(sample_data)[1], 0, 0.001) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)

D <- rdist(locc) # distance

alpha = 0.0587
beta  = 0.0520
delta = 0.0480


M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
S <- alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3))))   # spherical covariance matrix 
diag(S) = diag(S) + delta                     # nugget
Z = matrix(res0, ncol = 1)

B_estimate =solve(t(M) %*% (solve(S, tol = 5e-21) %*% M),tol= 5e-21) %*% t(M) %*% solve(S,tol = 5e-21)%*% Z
print(round(B_estimate,5))

##  BETA PARAMETERS FOR SPHERICAL OLS
#191.21560
#x1  -0.00013
#x2  -0.07771
#x3  -0.00182
#x4   3.01814
#x5   4.71503

##---------- Kriging  --------------------------------------------------------- 
res_krig <- (lm(log(median_house_value) ~ total_bedrooms+ 
                  median_income+housing_median_age + 
                  longitude+latitude,  data = sample_kriging_data))$residuals
#res_krig = matrix(res_krig, ncol = 1)

hist(res_krig)
# locations for selected data and kriging data
locc = cbind(sample_data$longitude, sample_data$latitude)
locc_k = cbind(sample_kriging_data$longitude, sample_kriging_data$latitude)


# distance
D <- rdist(locc)                  # distance b/t the selected data
d <- rdist(locc,  locc_k)         # distance b/t selected data and kriging data


k1     <- sample_kriging_data$total_bedrooms 
k2     <- sample_kriging_data$median_income  
k3     <- sample_kriging_data$housing_median_age 
k4     <- sample_kriging_data$longitude
k5     <- sample_kriging_data$latitude




# modify code below  
M = cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) 
m = t(cbind(rep(1,100), k1,k2,k3,k4,k5)) 

Z = matrix(res0, ncol = 1)       # data vector



# -------------------------SPHERICAL OLS RESULT
######## SPHERICAL OLS
alpha = 0.0587
beta  = 0.0520
delta = 0.0480

S <- alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3)))) # SPHERICAL COVARIANCE MATRIX 
diag(S) = diag(S)+delta          # add nugget 
k = alpha*(1-((3*d)/(2*beta))+((d**3)/(2*(beta**3))))
lambda = (solve(S, tol = 5e-21) - solve(S, tol = 5e-21) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-21) %*% M, tol = 5e-21) %*% t(M) %*% solve(S, tol = 5e-21) ) %*% k + solve(S, tol = 5e-21) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-21) %*% M, tol = 5e-21) %*% m 
dim(lambda)

krig.sph.ols = t(lambda) %*% Z 
dim(krig.sph.ols)




plot(res_krig, krig.sph.ols)
lines(res_krig,res_krig) #y=x line 

mse1=mean((krig.sph.ols-res_krig)^2)
mse1 #in log form
mae1=mean(abs(krig.sph.ols-res_krig))
mae1

#> mse1 #in log form
#[1] 0.02946355

#> mae1
#[1] 0.1173685


#################### END OF SPHERICAL OLS ON OVERALL  #########################

#################### START OF SPHERICAL WLS ON OVERALL  #########################

fit.sph.WLS=variofit(vario0, ini.cov.pars=c(0.15,0.05),
                     #weights = 'npairs',
                     cov.model = "spherical"
                     )
fit.sph.WLS
#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0000  0.1077  0.0390 
#Practical Range with cor=0.05 for asymptotic range: 0.03896988

#variofit: minimised weighted sum of squares = 139.7401

#NOTE: There is no nugget.
#Thus, this will lead to the covariance matrix being numerically singular


set.seed(222)
zz=rnorm(dim(sample_data)[1], 0, 0.001) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)

D <- rdist(locc) # distance

alpha = 0.1077
beta  = 0.0390
delta = 0.000


M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
S <- alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3))))   # spherical covariance matrix 

diag(S) = diag(S) + delta                     # nugget
Z = matrix(res0, ncol = 1)

B_estimate =solve(t(M) %*% (solve(S, tol = 5e-21) %*% M),tol= 5e-21) %*% t(M) %*% solve(S,tol = 5e-21)%*% Z
print(round(B_estimate,5))

##  BETA PARAMETERS FOR SPHERICAL WLS
#294.55132
#x1  -0.00001
#x2  -0.05928
#x3  -0.00231
#x4   3.50956
#x5   3.58601

##---------- Kriging  --------------------------------------------------------- 
res_krig <- (lm(log(median_house_value) ~ total_bedrooms+ 
                  median_income+housing_median_age + 
                  longitude+latitude,  data = sample_kriging_data))$residuals
#res_krig = matrix(res_krig, ncol = 1)

hist(res_krig)
# locations for selected data and kriging data
locc = cbind(sample_data$longitude, sample_data$latitude)
locc_k = cbind(sample_kriging_data$longitude, sample_kriging_data$latitude)


# distance
D <- rdist(locc)                  # distance b/t the selected data
d <- rdist(locc,  locc_k)         # distance b/t selected data and kriging data


k1     <- sample_kriging_data$total_bedrooms 
k2     <- sample_kriging_data$median_income  
k3     <- sample_kriging_data$housing_median_age 
k4     <- sample_kriging_data$longitude
k5     <- sample_kriging_data$latitude




# modify code below  
M = cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) 
m = t(cbind(rep(1,100), k1,k2,k3,k4,k5)) 

Z = matrix(res0, ncol = 1)       # data vector



# -------------------------SPHERICAL WLS RESULT
######## SPHERICAL WLS
alpha = 0.1077
beta  = 0.0390
delta = 0.000


S <- alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3)))) # SPHERICAL COVARIANCE MATRIX 
diag(S) = diag(S)+ 10^(-4)          # add nugget 
eigen(S,only.values = T)
k = alpha*(1-((3*d)/(2*beta))+((d**3)/(2*(beta**3))))
lambda = (solve(S, tol = 5e-25) - solve(S, tol = 5e-25) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-25) %*% M, tol = 5e-25) %*% t(M) %*% solve(S, tol = 5e-25) ) %*% k + solve(S, tol = 5e-25) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-25) %*% M, tol = 5e-25) %*% m 
dim(lambda)

krig.sph.wls = t(lambda) %*% Z 
dim(krig.sph.wls)




plot(res_krig, krig.sph.wls)
lines(res_krig,res_krig) #y=x line 

mse1=mean((krig.sph.wls-res_krig)^2)
mse1 #in log form
mae1=mean(abs(krig.sph.wls-res_krig))
mae1

#> mse1 #in log form
#[1] 0.04137684
#> mae1=mean(abs(krig.sph.wls-res_krig))
#> mae1
#[1] 0.123151

#################### END OF SPHERICAL WLS ON OVERALL  #########################


##################### END OF SPHERICAL WORK ON OVERALL DATA ###################



##########  START OF LOCAL STATIONARY WORK FOR SPHERICAL   #########################

# quick renaming 
data0  = sample_data
data.v = sample_kriging_data

range(sample_data$longitude)
data1=data0[data0$longitude < -122.445,] 
data2=data0[data0$longitude >= -122.445 & data0$longitude <  -122.300,]
data3=data0[data0$longitude >=  -122.300 & data0$longitude< -122.155,] 
data4=data0[data0$longitude >= -122.155,] 


par(mfrow=c(1,1),mai=c(0.5,0.5,0.5,0.5)) 
quilt.plot(data0$longitude, data0$latitude, res0) 
US(add=T) 
abline(v= c( -122.445, -122.300,-122.155),col="gray") # ablines at -99 and -97 longatude


text(-122.56,  37.6, label=nrow(data1), col="black") 
text(-122.4,   37.6, label=nrow(data2), col="black") 
text(-122.21,  37.6, label=nrow(data3), col="black") 
text(-122.12,  37.6, label=nrow(data4), col="black") 
print(
  cat(" Number of cases within each subgroup:", '\n',
      "data1 = ", nrow(data1) ,'\n',
      "data2 = ", nrow(data2) ,'\n',
      "data3 = ", nrow(data3) ,'\n',
      "data4 = ", nrow(data4) ,'\n'
  )
)

#nrow(data1);nrow(data2);nrow(data3);nrow(data4)

#################### START OF SPHERICAL OLS ON LOCAL STATIONARY  #########################

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
              max.dist = 1) 
plot(vario1)

#SPHERICAL 
plot(vario1, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.04,0.05), 
                 nugget = 0.005, 
                 max.dist = 1,
                 lwd = 3)

res_krig_data2 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data2))$residuals

vario2=variog(coords=cbind(data2$longitude, data2$latitude), 
              data=res_krig_data2, max.dist = 1,
              trend=trend.spatial(~data2$total_bedrooms  + 
                                    data2$median_income +
                                    data2$housing_median_age +
                                    data2$longitude+ 
                                    data2$latitude))
plot(vario2)
plot(vario2, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.10,0.05), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)


res_krig_data3 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data3))$residuals

vario3=variog(coords=cbind(data3$longitude, data3$latitude), 
              data=res_krig_data3, 
              max.dist = 0.45,
              trend=trend.spatial(~data3$total_bedrooms  + 
                                    data3$median_income +
                                    data3$housing_median_age +
                                    data3$longitude+ 
                                    data3$latitude
              ))
plot(vario3)
plot(vario3, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.10,0.10), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)

res_krig_data4 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data4))$residuals

vario4=variog(coords=cbind(data4$longitude, data4$latitude), 
              data=res_krig_data4, max.dist = 0.45,
              trend=trend.spatial(~data4$total_bedrooms  + 
                                    data4$median_income +
                                    data4$housing_median_age +
                                    data4$longitude+ 
                                    data4$latitude
              )) 
plot(vario4)
plot(vario4, main = "SPHERICAL" ) 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.05,0.05), 
                 nugget = 0.005, 
                 max.dist = 1,
                 lwd = 3)


#---plot-----------------------------------------------------

par(mfrow=c(2,2)) 
plot(vario1, main = "Data1")
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.04,0.05), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)

plot(vario2, main= "Data2") 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.10,0.05), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)
plot(vario3, main= "Data3")
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.10,0.10), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)
plot(vario4, main= "Data4") 
lines.variomodel(cov.model = "spherical", 
                 cov.pars = c(0.05,0.05), 
                 nugget = 0.005, 
                 max.dist = 1,
                 lwd = 3)
par(mfrow=c(1,1))

###------ SPHERICAL OLS  variofit() of the subgroups-------------------------------------

fit1=variofit(vario1, ini.cov.pars=c(0.04,0.1),
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



fit3=variofit(vario3, ini.cov.pars=c(0.1,0.1),
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

fit4=variofit(vario2, ini.cov.pars=c(0.05,0.05),
              weights = 'equal',
              cov.model = "spherical")
fit4

#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0682  0.0430  0.0621 
#Practical Range with cor=0.05 for asymptotic range: 0.06208332

#variofit: minimised sum of squares = 0.0223


#------- SPHERICAL OLS KRIGING FOR LOCAL STATIONARY

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
  
  lambda = (solve(S, tol = 5e-25) - solve(S, tol = 5e-25) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-25) %*% M, tol = 5e-25) %*% t(M) %*% solve(S, tol = 5e-25) ) %*% k + solve(S, tol = 5e-25) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-25) %*% M, tol = 5e-25) %*% m 
  
  krig=t(lambda) %*% lm(log(median_house_value) ~ total_bedrooms+ 
                          median_income+housing_median_age + 
                          longitude+latitude,  data = data)$residuals
  
  return(krig) 
} 


# ----   results 


## kriging with isotropic fit 
krig0 = krig(data0, data.v, c( 0.0587, 0.0520, 0.0480)) 

#--------------------------------------------------------change the range 

print(list(fit1, fit2, fit3, fit4))

## kriging over region 1 
data.v1 = data.v[data.v$longitude< -122.445,] 
krig2.1 = krig(data1,data.v1, c(0.0272,  0.0458, 0.0132)) 

## kriging over region 2 
data.v2 = data.v[data.v$longitude >= -122.445 & data.v$longitude <  -122.300,] 
krig2.2 = krig(data2,data.v2, c(0.0430,  0.0621, 0.0682  )) 

## kriging over region 3 
data.v3 = data.v[data.v$longitude >=  -122.300 & data.v$longitude< -122.155,] 
krig2.3 = krig(data3,data.v3, c(0.0619,  0.0324,0.0340 ))

## kriging over region 4 
data.v4 = data.v[data.v$longitude >= -122.155,] 
krig2.4 = krig(data4,data.v4, c(0.0430,  0.0682, 0.0621 )) 

# plotting the kriging 
par(mfrow=c(1,1)) 
plot(data0$longitude, data0$latitude, pch=20, main = "SPHERICAL OLS") 
US(add=T) 
abline(v=c(-122.445,-122.300,-122.155),col="gray") 

points(data.v1$longitude, data.v1$latitude, col=2, lwd = 3) 
points(data.v2$longitude, data.v2$latitude, col=3, lwd = 3) 
points(data.v3$longitude, data.v3$latitude, col=4, lwd = 3) 
points(data.v4$longitude, data.v4$latitude, col=5, lwd = 3)


# --------checking errors v true values----------------------------- 
par(mfrow=c(1,2))

quilt.plot(data.v$longitude, data.v$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                   median_income+housing_median_age + 
                                                   longitude+latitude,  data = data.v)$residuals - krig0,
           main= "overall") 
US(add=T) 
data.vv=rbind(data.v1, data.v2, data.v3, data.v4) 
krig.vv=c(krig2.1, krig2.2, krig2.3, krig2.4) 
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

fit1a=variofit(vario1, ini.cov.pars=c(0.04,0.1),
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



fit3a=variofit(vario3, ini.cov.pars=c(0.1,0.1),
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

fit4a=variofit(vario2, ini.cov.pars=c(0.05,0.05),
              weights = 'npairs',
              cov.model = "spherical")
fit4a

#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0680  0.0589  0.0877 
#Practical Range with cor=0.05 for asymptotic range: 0.08767602

#variofit: minimised weighted sum of squares = 70.0208


#------- SPHERICAL WLS KRIGING FOR LOCAL STATIONARY

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
  
  lambda = (solve(S, tol = 5e-25) - solve(S, tol = 5e-25) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-25) %*% M, tol = 5e-25) %*% t(M) %*% solve(S, tol = 5e-25) ) %*% k + solve(S, tol = 5e-25) %*% M %*% solve(t(M) %*% solve(S, tol = 5e-25) %*% M, tol = 5e-25) %*% m 
  
  krig=t(lambda) %*% lm(log(median_house_value) ~ total_bedrooms+ 
                          median_income+housing_median_age + 
                          longitude+latitude,  data = data)$residuals
  
  return(krig) 
} 


# ----   results 


## kriging with isotropic fit 
krig0a = krig(data0, data.v, c( 0.1077, 0.0390, 10^(-4)  )) 

#--------------------------------------------------------change the range 

print(list(fit1a, fit2a, fit3a, fit4a))

## kriging over region 1 
data.v1 = data.v[data.v$longitude< -122.445,] 
krig2.1a = krig(data1,data.v1, c(0.0321,  0.0558, 0.0133)) 

## kriging over region 2 
data.v2 = data.v[data.v$longitude >= -122.445 & data.v$longitude <  -122.300,] 
krig2.2a = krig(data2,data.v2, c(0.0591,  0.0873, 0.0678  )) 

## kriging over region 3 
data.v3 = data.v[data.v$longitude >=  -122.300 & data.v$longitude< -122.155,] 
krig2.3a = krig(data3,data.v3, c(0.0633,  0.0334,0.0342 ))

## kriging over region 4 
data.v4 = data.v[data.v$longitude >= -122.155,] 
krig2.4a = krig(data4,data.v4, c(0.0589,  0.0877, 0.0680 )) 

# plotting the kriging 
par(mfrow=c(1,1)) 
plot(data0$longitude, data0$latitude, pch=20, main = "SPHERICAL WLS") 
US(add=T) 
abline(v=c(-122.445,-122.300,-122.155),col="gray") 

points(data.v1$longitude, data.v1$latitude, col=2, lwd = 3) 
points(data.v2$longitude, data.v2$latitude, col=3, lwd = 3) 
points(data.v3$longitude, data.v3$latitude, col=4, lwd = 3) 
points(data.v4$longitude, data.v4$latitude, col=5, lwd = 3)


# --------checking errors v true values----------------------------- 
par(mfrow=c(1,2))

quilt.plot(data.v$longitude, data.v$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                   median_income+housing_median_age + 
                                                   longitude+latitude,  data = data.v)$residuals - krig0,
           main= "overall") 
US(add=T) 
data.vv=rbind(data.v1, data.v2, data.v3, data.v4) 
krig.vva=c(krig2.1a, krig2.2a, krig2.3a, krig2.4a) 
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