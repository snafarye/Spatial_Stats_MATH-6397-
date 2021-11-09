library(readxl)
load('C:/Users/nsara/Documents/Grad_School/MATH 6397/report/data1.RData')
dim(data1)
#library(fields)
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


# log transformation of mean_house_value? 
hist(sample_data$median_house_value)
hist(log(sample_data$median_house_value))


qqnorm(sample_data$median_house_value)
qqline(sample_data$median_house_value)

qqnorm(log(sample_data$median_house_value))
qqline(log(sample_data$median_house_value))

hist(res0);qqnorm(res0);qqline(res0)
# best to take the log of the median_house_value for better normality results
# later see that using the residual as the predicting value is the best option ? 


# variables to use for covariates variables ro help describe dependent variable ??
par(mfrow=c(2,3))
plot(sample_data$housing_median_age,      log(sample_data$median_house_value)) # no clear linear/non-linear trend 
plot(sample_data$total_bedrooms,          log(sample_data$median_house_value)) # some what of a pos linear trend, but not strong  
plot(sample_data$median_income,           log(sample_data$median_house_value)) # much more significant , pos trend

plot(sample_data$longitude,      log(sample_data$median_house_value)) # not really see a sig. relationship b/t house_median_value
plot(sample_data$latitude,       log(sample_data$median_house_value)) # and long./ lat 
par(mfrow=c(1,1))


# Use the same linear regression model you used in team report 1 (with the same covariates)
# parameters: total_bedrooms,   median_income,  housing_median_age

# linear regression model from report 1
fit0 = lm(log(median_house_value) ~ total_bedrooms+ 
           median_income+housing_median_age + 
            longitude+latitude,  data = sample_data)
summary(fit0)
res0 = fit0$residuals
summary(res0)

plot(res0)
quilt.plot(sample_data$longitude, sample_data$latitude, res0)
US(add=T) 
hist(res0);qqnorm(res0);qqline(res0)

par(mfrow=c(2,3))
plot(sample_data$housing_median_age,      res0) # no clear linear/non-linear trend 
plot(sample_data$total_bedrooms,          res0) # some what of a pos linear trend, but not strong  
plot(sample_data$median_income,           res0) # much more significant , pos trend

plot(sample_data$longitude,               res0) # not really see a sig. relationship b/t house_median_value
plot(sample_data$latitude,                res0) # and long./ lat 
par(mfrow=c(1,1))


# -----------  part 1-----------------------------------------------------------


## --- Exponential ===============================
loc = cbind(sample_data$longitude, sample_data$latitude)
y      <- log(sample_data$median_house_value) #value at loc
res0   =  fit0$residuals
x1     <- sample_data$total_bedrooms 
x2     <- sample_data$median_income  
x3     <- sample_data$housing_median_age 
x4     <- sample_data$longitude
x5     <- sample_data$latitude

#zz=rnorm(dim(sample_data)[1], 0, 0.001) 
#sample_data$longitudee=sample_data$longitude+zz 
#locc = cbind(sample_data$longitudee, sample_data$latitude)

vario0  = variog(coords=loc,
                     data=res0, 
                     trend=trend.spatial(~x1+x2+x3+x4+x5),
                max.dist = 0.6) 
plot(vario0)

par(mfrow=c(1,2))
# expodentail 
plot(vario0, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.1,0.05), 
                 nugget = 0.01, 
                 max.dist = 0.6,
                 lwd = 3)

# spherical 
plot(vario0 , main= "SPH") 
lines.variomodel(cov.model = "sph", 
                 cov.pars = c(0.1,0.2), 
                 nugget = 0.01, 
                 max.dist = 0.6,  lwd = 3)
par(mfrow=c(1,1))



#-----same exact variogram plot -------------------
plot(vario0, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.1,0.05), 
                 nugget = 0.01, 
                 max.dist = 0.6,
                 lwd = 3)

# ------------------OLS
fit.expo.OLS=variofit(vario0, ini.cov.pars=c(0.10,0.05),
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
  
# ----------------WLS

fit.expo.WLS=variofit(vario0, ini.cov.pars=c(0.10,0.05),
                      cov.model = "exponential") 
fit.expo.WLS
#variofit: model parameters estimated by WLS (weighted least squares):
##  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.0332  0.0746  0.0188 
#Practical Range with cor=0.05 for asymptotic range: 0.05629127
#variofit: minimised weighted sum of squares = 137.8366

alpha = 0.0746
beta  = 0.0188
delta = 0.0332
  
  
  M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
  S <- alpha*exp(-D/beta)                       # covareiance matrix 
  diag(S) = diag(S) + delta                     # nugget
  Z = matrix(res0, ncol = 1)
  
  B_estimate2 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
  print(round(B_estimate2,5))
  
  
  50.01568
  x1 -0.00012
  x2 -0.07073
  x3 -0.00217
  x4  0.46468
  x5  0.19241
  



#--------------- REML
#likelihood method

fit.exp.reml=likfit(coords=locc, 
                 data=res0, 
                 trend=trend.spatial(~x1+x2+x3+x4+x5), # mean structure
                 lik.method = "REML",
                 cov.model = "exponential",
                 ini.cov.pars=c(0.10,0.05)) 
fit.exp.reml

#likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
#"63.7370" "-0.0001" "-0.0749" "-0.0019" " 0.6377" " 0.3898" " 0.0312" " 0.1026" " 0.0796" 
#Practical Range with cor=0.05 for asymptotic range: 0.2384918
#likfit: maximised log-likelihood = 124.5

plot(vario0, main  = "Expo")
lines(fit.expo.OLS, lwd = 2, col = "green")
lines(fit.expo.WLS, lwd = 2, col= "red")
lines(fit.exp.reml, lwd = 2, col = "blue")
legend(0.3,0.02, legend=c("OLS","WLS", "REML"),
       lwd=c(2,2,2),col = c("green","red","blue"), cex=0.7)


#==================================================================
#==================================================================
#==================================================================

## --- spherical ===============================

# reminder 
vario0 = variog(coords=loc,
                data=res0, 
                trend=trend.spatial(~x1+x2+x3+x4+x5), 
                max.dist = 0.6) 
plot(vario0) 

# spherical 
plot(vario0 , main= "SPH") 
lines.variomodel(cov.model = "sph", 
                 cov.pars = c(0.10,0.2), 
                 nugget = 0.02, 
                 max.dist = 0.6,  lwd = 3)


# ------------------OLS

fit.sph.OLS=variofit(vario0, ini.cov.pars=c(0.10,0.2),
                      weights = 'equal',
                     cov.model = "spherical") ## try many other models 
fit.sph.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0480  0.0587  0.0520 
#Practical Range with cor=0.05 for asymptotic range: 0.05196646
#variofit: minimised sum of squares = 0.0016


alpha = 0.0587
beta  = 0.0520
delta = 0.0480

  
  M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5)                 # design matrix
  S = alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3))))         # covareiance matrix 
  diag(S) = diag(S) + delta                                     # nugget
  Z = matrix(res0, ncol = 1)
  
  B_estimate4 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
  print(round(B_estimate4,5))
  

    80.22697
  x1 -0.00013
  x2 -0.07776
  x3 -0.00182
  x4  2.23916
  x5  5.13734
  
# ----------------WLS

fit.sph.WLS=variofit(vario0, ini.cov.pars=c(0.10,0.2),
                      weights = 'npairs',
                     cov.model = "spherical") ## try many other models 
fit.sph.WLS

#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0481  0.0596  0.0531 
#Practical Range with cor=0.05 for asymptotic range: 0.05312656
#variofit: minimised weighted sum of squares = 136.576



alpha = 0.0596
beta  = 0.0531
delta = 0.0481


M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5)                 # design matrix
S = alpha*(1-((3*D)/(2*beta)) + ((D^3)/(2*(beta^3))))  # covareiance matrix 
diag(S) = diag(S) + delta                                     # nugget
Z = matrix(res0, ncol = 1)

B_estimate5 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z

print(round(B_estimate5,5))

84.12336
x1 -0.00013
x2 -0.07765
x3 -0.00183
x4  2.25262
x5  5.07744

#--------------- REML

fit.sph.reml=likfit(coords=locc, 
                    data=res0, 
                    trend=trend.spatial(~x1+x2+x3+x4+x5), # mean structure
                    lik.method = "REML",
                    cov.model = "spherical",
                    ini.cov.pars=c(0.10,0.2)) 
fit.sph.reml

#likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
#"59.1890" "-0.0001" "-0.0754" "-0.0018" " 0.5860" " 0.3426" " 0.0318" " 0.1005" " 0.1321" 
#Practical Range with cor=0.05 for asymptotic range: 0.132113
#likfit: maximised log-likelihood = 123.8



plot(vario0, main = "Spherical")
lines(fit.sph.OLS,  lwd = 2, col = "green")
lines(fit.sph.WLS,  lwd = 2, col= "red")
lines(fit.sph.reml, lwd = 2, col = "blue")
legend(0.3,0.02, legend=c("OLS","WLS", "REML"),
       lwd=c(2,2,2),col = c("green","red","blue"), cex=0.7)






