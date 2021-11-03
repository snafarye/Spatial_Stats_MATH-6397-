library(readxl)
load('C:/Users/nsara/Documents/Grad_School/MATH 6397/report/data1.RData')
dim(data1)
library(fields)
library(gstat) 
library(sp) 
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
fit0 = lm(log(median_house_value) ~longitude+latitude+
           total_bedrooms+ 
           median_income+housing_median_age, 
          data = sample_data)
summary(fit0)
res0 = fit0$residuals
plot(res0)
quilt.plot(sample_data$longitude, sample_data$latitude, res0)

hist(res0);qqnorm(res0);qqline(res0)

par(mfrow=c(2,3))
plot(sample_data$housing_median_age,      res0) # no clear linear/non-linear trend 
plot(sample_data$total_bedrooms,          res0) # some what of a pos linear trend, but not strong  
plot(sample_data$median_income,           res0) # much more significant , pos trend

plot(sample_data$longitude,               res0) # not really see a sig. relationship b/t house_median_value
plot(sample_data$latitude,                res0) # and long./ lat 
par(mfrow=c(1,1))

loc = cbind(sample_data$longitude, sample_data$latitude)

plot(vgram(loc,log(sample_data$median_house_value)),lon.lat=TRUE) #Signs of spacial varying mean
plot(vgram(loc,log(sample_data$median_house_value),dmax = 0.6), lon.lat=TRUE)

plot(vgram(loc,sample_data$median_house_value),lon.lat=TRUE) #still have Signs of spacial varying mean
plot(vgram(loc,res0),lon.lat=TRUE)  # use residual to get a better variogram plot
plot(vgram(loc,res0, dmax = 0.5),lon.lat=TRUE)  # set a dmax for the plot

# -----------  part 1-----------------------------------------------------------


## --- Exponential ===============================
loc = cbind(sample_data$longitude, sample_data$latitude)
y      <- log(sample_data$median_house_value) #value at loc
y_res  <- res0
x1     <- sample_data$total_bedrooms 
x2     <- sample_data$median_income  
x3     <- sample_data$housing_median_age 
x4     <- sample_data$longitude
x5     <- sample_data$latitude

vario1 = variog(coords=loc,
              data=y, 
              trend=trend.spatial(~x1+x2+x3+x4+x5) ) 
plot(vario1) 

#vario.y_res = variog(coords=loc,
#                     data=y_res, 
#                     trend=trend.spatial(~x1+x2+x3+x4+x5)) 
#plot(vario.y_res) 
#-----same exact variogram plot -------------------

# ------------------OLS

fit.expo.OLS=variofit(vario1, ini.cov.pars=c(0.15,0.2),
                      weights = 'equal',
                      cov.model = "exponential") ## try many other models 
fit.expo.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.0479  0.0542  0.0186 
#Practical Range with cor=0.05 for asymptotic range: 0.05568734
#variofit: minimised sum of squares = 0.0017


D <- rdist.earth(locc,miles = T) # distance

# take cov para estimator from OLS variofit() above
alpha = 0.0542
beta  = 0.0186
delta = 0.0479


OLS.para= function(par){
  print(par)
  alpha = exp(par[1]) 
  beta  = exp(par[2])
  delta = exp(par[3])
  nu= exp(par[4])/(1+exp(par[4]))*5 
    
  
  M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
  S <- alpha*exp(-D/beta)                       # covareiance matrix 
  diag(S) = diag(S) + delta                     # nugget
  Z = matrix(y, ncol = 1)
  
  B_estimate = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z

  print(B_estimate)
  return(B_estimate)
  
  }
ini = c(0.0542,0.0186,0.0479)
fit.OLS <- optim(ini, OLS.para
                 #, control= list(maxit= 100000000)
                 )
[,1]
-1.080342e+02
x1  1.808418e-04
x2  1.402975e-01
x3  1.927285e-03
x4 -1.197525e+00
x5 -7.083061e-01

# ----------------WLS

fit.expo.WLS=variofit(vario1, ini.cov.pars=c(0.15,0.2),
                      weights = 'npairs',
                      cov.model = "exponential") ## try many other models 
fit.expo.WLS
#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.1058  0.0000  0.1792 
#Practical Range with cor=0.05 for asymptotic range: 0.5367151
#variofit: minimised weighted sum of squares = 139.1105

D <- rdist.earth(locc,miles = T) # distance
# D <- rdist(locc,miles =T)
alpha = 0.0000
beta  = 0.1792
delta = 0.1058


WLS.para= function(par){
  print(par)
  alpha = exp(par[1]) 
  beta  = exp(par[2])
  delta = exp(par[3])
  nu= exp(par[4])/(1+exp(par[4]))*5 
  
  
  M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
  S <- alpha*exp(-D/beta)                       # covareiance matrix 
  diag(S) = diag(S) + delta                     # nugget
  Z = matrix(y, ncol = 1)
  
  B_estimate = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
  
  print(B_estimate)
  return(B_estimate)
  
}
ini = c(alpha, beta, delta)
fit.WLS <- optim(ini, 
                 WLS.para
                 #, control= list(maxit= 100000000)
)
   -1.087199e+02
x1  1.786055e-04
x2  1.397721e-01
x3  1.855820e-03
x4 -1.202295e+00
x5 -7.054010e-01

#--------------- REML
#likelihood method
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)


fit.exp.reml=likfit(coords=locc, 
                 data=y, 
                 trend=trend.spatial(~x1+x2+x3+x4+x5), # mean structure
                 lik.method = "REML",
                 cov.model = "exponential",
                 ini.cov.pars=c(0.15,0.2)) 
fit.exp.reml

#likfit: estimated model parameters:
#  beta0       beta1       beta2       beta3       beta4       beta5       tausq     sigmasq         phi 
#"-127.3415" "   0.0002" "   0.1413" "   0.0020" "  -1.3370" "  -0.6487" "   0.0762" "   0.0375" "   0.0905" 
#Practical Range with cor=0.05 for asymptotic range: 0.2711601
#likfit: maximised log-likelihood = -325.1

#------ with data = y_res, notice only the parameters change

#likfit: estimated model parameters:
#  beta0      beta1      beta2      beta3      beta4      beta5      tausq    sigmasq        phi 
#"-22.7239" "  0.0000" " -0.0119" " -0.0010" " -0.1681" "  0.0602" "  0.0762" "  0.0375" "  0.0905" 
#Practical Range with cor=0.05 for asymptotic range: 0.2711601
#likfit: maximised log-likelihood = -325.1



#==================================================================

# changing the aplha and beta values to 0.30 and 0.2
#likelihood method
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)


fit.exp.reml=likfit(coords=locc, 
                 data=y, 
                 trend=trend.spatial(~x1+x2+x3+x4+x5), # mean structure
                 lik.method = "REML",
                 cov.model = "exponential",
                 ini.cov.pars=c(0.30,0.2)) 
fit.exp.reml
#likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     tausq   sigmasq       phi 
#"11.4908" " 0.0002" " 0.1446" " 0.0037" " 0.0871" " 0.0593" " 0.2020" 
#Practical Range with cor=0.05 for asymptotic range: 0.6050156
#likfit: maximised log-likelihood = -401.1



## --- spherical ===============================

# reminder 
vario1 = variog(coords=loc,
                data=y, 
                trend=trend.spatial(~x1+x2+x3+x4+x5) ) 
plot(vario1) 

# ------------------OLS

fit.sph.OLS=variofit(vario2, ini.cov.pars=c(0.15,0.2),
                      weights = 'equal',
                     cov.model = "spherical") ## try many other models 
fit.sph.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0539  0.0761  0.1027 
#Practical Range with cor=0.05 for asymptotic range: 0.1026624
#variofit: minimised sum of squares = 0.0019


# ----------------WLS

fit.sph.WLS=variofit(vario2, ini.cov.pars=c(0.15,0.2),
                      weights = 'npairs',
                     cov.model = "spherical") ## try many other models 
fit.sph.WLS
#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0540  0.0846  0.1156 
#Practical Range with cor=0.05 for asymptotic range: 0.1155594
#variofit: minimised weighted sum of squares = 64.0742

#--------------- REML

fit.sph.reml=likfit(coords=locc, 
                    data=y, 
                    trend=trend.spatial(~x1+x2+x3+x4+x5), # mean structure
                    lik.method = "REML",
                    cov.model = "spherical",
                    ini.cov.pars=c(0.15,0.2)) 
fit.sph.reml
#likfit: estimated model parameters:
#  beta0       beta1       beta2       beta3       beta4       beta5       tausq     sigmasq         phi 
#"-129.9206" "   0.0002" "   0.1412" "   0.0020" "  -1.3486" "  -0.6180" "   0.0779" "   0.0664" "   0.3501" 
#Practical Range with cor=0.05 for asymptotic range: 0.3501326
#likfit: maximised log-likelihood = -326.6


#----------- with data =  y_res, notice the paramerter are diff, but the max log-lik is the same
#likfit: estimated model parameters:
#  beta0      beta1      beta2      beta3      beta4      beta5      tausq    sigmasq        phi 
#"-25.3030" "  0.0000" " -0.0120" " -0.0009" " -0.1797" "  0.0909" "  0.0779" "  0.0664" "  0.3501" 
#Practical Range with cor=0.05 for asymptotic range: 0.3501326
#likfit: maximised log-likelihood = -326.6


















