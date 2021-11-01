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
# best to take the log of the median_house_value for better normality results


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

loc = cbind(sample_data$longitude, sample_data$latitude)

plot(vgram(loc,log(sample_data$median_house_value)),lon.lat=TRUE) #Signs of spacial varying mean
plot(vgram(loc,log(sample_data$median_house_value),dmax = 0.6), lon.lat=TRUE)

plot(vgram(loc,sample_data$median_house_value),lon.lat=TRUE) #still have Signs of spacial varying mean
plot(vgram(loc,res0),lon.lat=TRUE)  # use residual to get a better variogram plot
plot(vgram(loc,res0, dmax = 0.5),lon.lat=TRUE)  # set a dmax for the plot

# -----------  part 1-----------------------------------------------------------
# Fill in the following table with estimates using the estimation method 
# mentioned below. (p: the number of covariates.)

loc = cbind(sample_data$longitude, sample_data$latitude)
y <- log(sample_data$median_house_value) #value at loc
x1  <- sample_data$total_bedrooms #x_total_bedrooms
x2  <- sample_data$median_income  # x_median_income
x3  <- sample_data$housing_median_age 
x4  <- sample_data$longitude
x5  <- sample_data$longitude


# field
plot(vgram(loc,y, lon.lat = FALSE))
plot(vgram(loc,y, lon.lat = TRUE, dmax = 30))

# geoR
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
fit = likfit(coords= cbind(sample_data$longitudee,sample_data$latitude),
            data = log(sample_data$median_house_value), 
            #trend=trend.spatial(~sample_data$total_bedrooms), # mean structure
            lik.method = "ML",
            ini.cov.pars=c(1,1)) 
fit
#likfit: estimated model parameters:
#  beta     tausq   sigmasq       phi 
#"12.4667" " 0.1436" " 0.1255" " 0.1713" 
#Practical Range with cor=0.05 for asymptotic range: 0.5133054
#likfit: maximised log-likelihood = -811.3
##--------------------------------------------------------------------------

fit.geo =variofit(variog(coords = loc, data=y,
                         trend = '1st'),
                  weights = "equal",
                  fix.kappa=FALSE, 
                  ini.cov.pars=c(0.25,0.2)) 
fit.geo
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: matern
#parameter estimates:
#  tausq sigmasq     phi   kappa 
#0.0425  0.1284  0.0024  1.0000 
#Practical Range with cor=0.05 for asymptotic range: 0.009516808
#variofit: minimised sum of squares = 0.0223
#------------------------------------------------------

vario0 <- variog(coords = loc, data=res0,
                trend = '1st',
               )
plot(vario0)
## variogram  (original data, not the res) 
vario1=variog(coords=loc, 
              data=y) 
plot(vario1) 

vario2=variog(coords=loc,
              data=y, 
              trend=trend.spatial(~x1+x2+x3) 
              ) 
plot(vario2) 
# least squares 
fit.expo=variofit(vario2, ini.cov.pars=c(0.15,0.2)) ## try many other models 
fit.expo


#------------------------------------------------------
lines.variomodel(cov.model="exp", cov.pars=c(0.3,.2), nug=0, max.dist=1)
ml <- likfit(coords = locc,
             data = y,
             ini = c(0.15,0.2))
ml


reml <- likfit(locc, data = y,ini = c(0.15,0.2), fix.nugget = T, method = "RML")
ols <- variofit(locc,data = y, ini = c(0.15,0.2), fix.nugget = TRUE, weights="equal")
ols
wls <- variofit(locc, data = y,ini = c(0.15,0.2), fix.nugget = T)


#gstat


coordinates(sample_data) = ~ longitude + latitude
vg = variogram(res0~loc, data = sample_data)
plot(vg)
model = vgm(psill = 0.25, 
            model="Exp", 
            nugget=0.01, 
            range=0.1)
plot(vg,model)
fit = fit.variogram(vg,model)
fit
plot(vg,fit)


coordinates(sample_data) = ~ longitude + latitude
vg.og = variogram(y~loc, data = sample_data)
plot(vg.og)
model = vgm(psill = 0.25, 
            model="Exp", 
            nugget=0.01, 
            range=0.15)
plot(vg.og,model)
fit = fit.variogram(vg.og,model)
fit
plot(vg.og,fit)



## --- Exponential ===============================
loc = cbind(sample_data$longitude, sample_data$latitude)
y <- log(sample_data$median_house_value) #value at loc

x1  <- sample_data$total_bedrooms #x_total_bedrooms
x2  <- sample_data$median_income  # x_median_income
x3  <- sample_data$housing_median_age 
x4  <- sample_data$longitude
x5  <- sample_data$longitude

vario2=variog(coords=loc,
              data=y, 
              trend=trend.spatial(~x1+x2+x3) 
) 
plot(vario2) 



# ------------------OLS

fit.expo.OLS=variofit(vario2, ini.cov.pars=c(0.15,0.2),
                      weights = 'equal',
                      cov.model = "exponential") ## try many other models 
fit.expo.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.0528  0.0775  0.0431 
#Practical Range with cor=0.05 for asymptotic range: 0.1291828
#variofit: minimised sum of squares = 0.002



# ----------------WLS

fit.expo.WLS=variofit(vario2, ini.cov.pars=c(0.15,0.2),
                      weights = 'npairs',
                      cov.model = "exponential") ## try many other models 
fit.expo.WLS
#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.0342  0.1056  0.0450 
#Practical Range with cor=0.05 for asymptotic range: 0.1348015
#variofit: minimised weighted sum of squares = 56.0108


#--------------- REML
#likelihood method
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)

fit2.exp.ML=likfit(coords=locc, 
            data=y, 
            trend=trend.spatial(~x1+x2+x3), # mean structure
            lik.method = "ML",
            cov.model = "exponential",
            ini.cov.pars=c(0.15,0.2))  
fit2.exp.ML
#likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     tausq   sigmasq       phi 
#"11.0976" " 0.0003" " 0.1750" " 0.0100" " 0.1254" " 0.0088" " 0.0082" 
#Practical Range with cor=0.05 for asymptotic range: 0.02471925

#likfit: maximised log-likelihood = -620.9



fit2.exp.reml=likfit(coords=locc, 
                 data=y, 
                 trend=trend.spatial(~x1+x2+x3), # mean structure
                 lik.method = "REML",
                 cov.model = "exponential",
                 ini.cov.pars=c(0.15,0.2)) 
fit2.exp.reml

#likfit: estimated model parameters:
#  beta0      beta1      beta2      beta3      tausq    sigmasq 
#" 11.0971" "  0.0003" "  0.1750" "  0.0100" "  0.1341" "  0.7424" 
#phi 
#"615.8160" 
#Practical Range with cor=0.05 for asymptotic range: 1844.82

#likfit: maximised log-likelihood = -621.6

#==================================================================

# changing the aplha and beta values to 0.30 and 0.2
#likelihood method
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)

fit.exp.ML=likfit(coords=locc, 
            data=y, 
            trend=trend.spatial(~x1+x2+x3), # mean structure
            lik.method = "ML",
            cov.model = "exponential",
            ini.cov.pars=c(0.30,0.2))  
fit.exp.ML
#does not run with the X4 and X5
##likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     tausq   sigmasq       phi 
#"11.4863" " 0.0002" " 0.1447" " 0.0037" " 0.0867" " 0.0522" " 0.1719" 
#Practical Range with cor=0.05 for asymptotic range: 0.5149146
#likfit: maximised log-likelihood = -402.8

fit.exp.reml=likfit(coords=locc, 
                 data=y, 
                 trend=trend.spatial(~x1+x2+x3), # mean structure
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

vario2=variog(coords=loc,
              data=y, 
              trend=trend.spatial(~x1+x2+x3) 
) 
plot(vario2) 



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
                    trend=trend.spatial(~x1+x2+x3), # mean structure
                    lik.method = "REML",
                    cov.model = "spherical",
                    ini.cov.pars=c(0.15,0.2)) 
fit.sph.reml
#likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     tausq   sigmasq       phi 
#"11.4706" " 0.0002" " 0.1446" " 0.0036" " 0.0877" " 0.0603" " 0.3603" 
#Practical Range with cor=0.05 for asymptotic range: 0.3603041
#likfit: maximised log-likelihood = -400.6






















