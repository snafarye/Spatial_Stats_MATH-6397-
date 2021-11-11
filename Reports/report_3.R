library(readxl)
load('C:/Users/nsara/Documents/Grad_School/MATH 6397/report/data1.RData')
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

#------------------- variogram----------------
vario0  = variog(coords=loc,
                 data=res0, 
                 trend=trend.spatial(~x1+x2+x3+x4+x5),
                 max.dist = 0.6) 
plot(vario0) 

# expodentail 
plot(vario0, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.1,0.05), 
                 nugget = 0.01, 
                 max.dist = 0.6,
                 lwd = 3)

# ------------------ OLS ------------------------------------------------

fit.expo.OLS=variofit(vario0, ini.cov.pars=c(0.15,0.05),
                      weights = 'equal',
                      cov.model = "exponential")
fit.expo.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.0474  0.0594  0.0210 
#Practical Range with cor=0.05 for asymptotic range: 0.06277277
#variofit: minimised sum of squares = 0.0017


set.seed(222)
zz=rnorm(dim(sample_data)[1], 0, 0.001) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)

D <- rdist(locc) # distance

alpha = 0.0594
beta  = 0.0210
delta = 0.0474

M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
S <- alpha*exp(-D/beta)                       # covareiance matrix 
diag(S) = diag(S) + delta                     # nugget
Z = matrix(res0, ncol = 1)

B_estimate1 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
print(round(B_estimate1,5))


48.62891
x1 -0.00012
x2 -0.06489
x3 -0.00212
x4  0.45398
x5  0.19365

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


# modify code below  
M = cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) 
m = t(cbind(rep(1,100), k1,k2,k3,k4,k5)) 

Z = matrix(res0, ncol = 1)       # data vector

# -------------------------(my results)
######## exponential OLS
alpha = 0.0594
beta  = 0.0210
delta = 0.0474

S = alpha * exp(-D/beta)         # covariance matrix 
diag(S) = diag(S)+delta          # add nugget 
k = alpha*exp(-d/beta)
lambda = (solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
dim(lambda)

krig.exp.ols = t(lambda) %*% Z 
dim(krig.exp.ols)




plot(res_krig, krig.exp.ols)
lines(res_krig,res_krig) #y=x line 

mse1=mean((krig.exp.ols-res_krig)^2)
mse1 #in log form
mae1=mean(abs(krig.exp.ols-res_krig))
mae1



> mse1=mean((krig.exp.ols-res_krig)^2)
> mse1 #in log form
[1] 0.03133114
> mae1=mean(abs(krig.exp.ols-res_krig))
> mae1
[1] 0.1226097



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

zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)


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

plot(vario1, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.08,0.02), 
                 nugget = 0.005, 
                 max.dist = 0.6,
                 lwd = 3)

fit1=variofit(vario1, ini.cov.pars=c(0.08,0.02),
              weights = 'equal',
              cov.model = "exponential")
fit1


res_krig_data2 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data2))$residuals

vario2=variog(coords=cbind(data2$longitude, data2$latitude), 
              data=res_krig_data2, max.dist = 0.6,
              trend=trend.spatial(~data2$total_bedrooms  + 
                                    data2$median_income +
                                    data2$housing_median_age +
                                    data2$longitude+ 
                                    data2$latitude))
plot(vario2)
plot(vario2, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.10,0.05), 
                 nugget = 0.01, 
                 max.dist = 1,
                 lwd = 3)

fit2=variofit(vario2, ini.cov.pars=c(0.1,0.05),
              weights = 'equal',
              cov.model = "exponential")
fit2

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
plot(vario3, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.08,0.05), 
                 nugget = 0.005, 
                 max.dist = 0.39,
                 lwd = 3)


fit3=variofit(vario3, ini.cov.pars=c(0.08,0.05),
              weights = 'equal',
              cov.model = "exponential")
fit3 
print(fit1;fit2;fit3)

#------estimate paramerters-------------------------------------------------
b_estimates = function(data, par){
  n=dim(data)[1] 
  alpha=par[1] 
  beta=par[2] 
  delta=par[3] 
  M=cbind(rep(1,n), data$total_bedrooms,data$median_income,data$housing_median_age,data$longitude,data$latitude) 
  zz=rnorm(dim(data)[1], 0, 0.001) 
  data$longitudee=data$longitude+zz 
  locc = cbind(data$longitudee, data$latitude)
  
  
  D=rdist(locc) 
  S=alpha*exp(-D/beta) 
  diag(S)=diag(S)+delta 
  
  
  b_estimates=solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% matrix(lm(log(median_house_value) ~ total_bedrooms+ 
                                                                                     median_income+housing_median_age + 
                                                                                     longitude+latitude,  data = data)$residuals, ncol = 1)
  
  return((round(b_estimates,5)))
  
}

b_estimates(data1, c(0.0521,0.1048, 0.0341 ))
b_estimates(data2, c(0.0738  ,0.0172 , 0.0349 ))
b_estimates(data3, c(0.0255  ,0.0145 , 0.0555 ))

#------- kriging the local stay 

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
  S=alpha*exp(-D/beta) 
  diag(S)=diag(S)+delta 
  
 
  
  k=alpha*exp(-d/beta) 
  
  lambda=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m
  krig=t(lambda) %*% lm(log(median_house_value) ~ total_bedrooms+ 
                          median_income+housing_median_age + 
                          longitude+latitude,  data = data)$residuals
  
  return(krig) 
} 


# ----   results 


## kriging with isotropic fit 
krig0 = krig(data0, data.v, c( 0.0542, 0.0186, 0.0479)) 



#--------------------------------------------------------change the range 

## kriging over region 1 
data.v1 = data.v[data.v$longitude > -122.204,] 
#krig2.1 = krig(data1,data.v1, c(0.0274,  0.0146, 0.0128)) 

## kriging over region 2 
data.v2 = data.v[data.v$longitude <= -122.204 & data.v$longitude>  -122.397,] 
#krig2.2 = krig(data2,data.v2, c(0.0435,  0.0189, 0.0674  )) 

## kriging over region 3 
data.v3 = data.v[data.v$longitude <=  -122.397 & data.v$longitude > -122.59,] 
#krig2.3 = krig(data3,data.v3, c(0.0621,  0.0123,0.0338 ))



par(mfrow=c(1,1)) 
plot(data0$longitude, data0$latitude, pch=20, main = "EXPO OLS") 
US(add=T) 
abline(v=c(-122.397, -122.204),col="gray") 

points(data.v1$longitude, data.v1$latitude, col=2, lwd = 3) 
points(data.v2$longitude, data.v2$latitude, col=3, lwd = 3) 
points(data.v3$longitude, data.v3$latitude, col=4, lwd = 3) 


# --------checking errors v true values----------------------------- 

par(mfrow=c(1,2))

# data.v --> overall stationarity 
quilt.plot(data.v$longitude, data.v$latitude, lm(log(median_house_value) ~ total_bedrooms+ 
                                                   median_income+housing_median_age + 
                                                   longitude+latitude,  data = data.v)$residuals - krig0,
           main= "overall") 
US(add=T) 

# data.vv --> local stationarity 
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





