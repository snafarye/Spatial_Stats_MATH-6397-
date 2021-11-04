library(readxl)
load("C:/Users/Samue/Desktop/books for school/FALL 2021/MATH 6397/data1.RData")
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



#LINEAR REGRESSION
fit0 = lm(log(median_house_value) ~longitude+latitude+
            total_bedrooms+ 
            median_income+housing_median_age, 
          data = sample_data)
summary(fit0)
res0 = fit0$residuals
plot(res0)

#QUILT PLOT
quilt.plot(sample_data$longitude, sample_data$latitude, res0)


#LOCATION STRUCTURING
loc = cbind(sample_data$longitude, sample_data$latitude)
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)

#COVARIATES
y      <- log(sample_data$median_house_value) #value at loc
y_res  <- res0
x1     <- sample_data$total_bedrooms 
x2     <- sample_data$median_income  
x3     <- sample_data$housing_median_age 
x4     <- sample_data$longitude
x5     <- sample_data$latitude

#VARIOGRAM
vario1 = variog(coords=loc,
                data=y_res, 
                trend=trend.spatial(~x1+x2+x3+x4+x5) ) 
plot(vario1) 



#ESTIMATION METHODS - SPHERICAL OLS
ols <- variofit(vario1,
                weights = "equal",
                #fix.kappa=FALSE, 
                cov.model = "spherical",
                ini.cov.pars=c(0.15,0.2)) 


ols

#RESULT FOR OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: spherical
#parameter estimates:
 # tausq sigmasq     phi 
#0.0480  0.0542  0.0582 
#Practical Range with cor=0.05 for asymptotic range: 0.05819513

#variofit: minimised sum of squares = 0.0017



#COVARIANCE PARAMETER ESTIMATES FROM OLS variofit()
D <- rdist.earth(locc,miles = T) # distance

alpha = 0.0542
beta  = 0.0582
delta = 0.0480

M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
#S <- alpha*exp(-D/beta)       # covareiance matrix 
S = 00542*(1-((3*D)/(2*0.0582))+((D**3)/(2*(0.0582**3)))) 
diag(S) = diag(S) + 0.0480                  # nugget
z=matrix(res0, ncol = 1)

B_estimate=solve(t(M) %*% (solve(S, tol = 5e-21) %*% M),tol= 5e-21) %*% t(M) %*% solve(S,tol = 5e-21)%*% z
print(B_estimate)
##RESULT OF B_Estimate:
#print(B_estimate)
#[,1]
#-1.615263e+03
#x1  4.364464e-04
#x2  3.793987e-02
#x3  1.225837e-02
#x4 -9.513237e-01
#x5  4.085006e+01


#ESTIMATION METHODS - SPHERICAL WLS
wls <- variofit(vario1,
                # weights = "equal",
                #fix.kappa=FALSE, 
                cov.model = "spherical",
                ini.cov.pars=c(0.15,0.2)) 
wls

#RESULT FOR WLS
#variofit: model parameters estimated by WLS (weighted least squares):
#  covariance model is: spherical
#parameter estimates:
#  tausq sigmasq     phi 
#0.0480  0.0593  0.0663 
#Practical Range with cor=0.05 for asymptotic range: 0.06629026

#variofit: minimised weighted sum of squares = 116.589

#PARAMETER ESTIMATION FUNCTION FOR WLS

D <- rdist.earth(locc,miles = T) # distance

#COVARIANCE PARAMETER ESTIMATES FROM WLS variofit()
alpha = 0.0593
beta  = 0.0663
delta = 0.0480

M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
#S <- alpha*exp(-D/beta)       # covareiance matrix 
S = 00593*(1-((3*D)/(2*0.0663))+((D**3)/(2*(0.0663**3)))) 
diag(S) = diag(S) + 0.0480                  # nugget
z=matrix(res0, ncol = 1)

B_estimate2=solve(t(M) %*% (solve(S, tol = 5e-21) %*% M),tol= 5e-21) %*% t(M) %*% solve(S,tol = 5e-21)%*% z
print(B_estimate2)
##RESULT FOR B ESTIMATE 2
#[,1]
#-1.611375e+03
#x1  4.165170e-04
#x2  3.469347e-02
#x3  1.186461e-02
#x4 -9.306397e-01
#x5  4.075791e+01

#ESTIMATION METHODS - SPHERICAL REML

reml <- likfit(coords= locc,
               data = y_res, 
               trend=trend.spatial(~x1+x2+x3+x4+x5), 
               lik.method = "REML",
               cov.model = "spherical",
               ini.cov.pars=c(0.15,0.2)) 

#RESULT FOR REML
#> reml
#likfit: estimated model parameters:
#  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
#"-9.9918" " 0.0000" "-0.0137" "-0.0012" "-0.0663" " 0.0533" " 0.0761" " 0.0477" " 0.2225" 
#Practical Range with cor=0.05 for asymptotic range: 0.2224812

#likfit: maximised log-likelihood = -318.5


#COVARIANCE SPHERICAL FUNCTION
S = alpha*(1-((3*D)/(2*beta))+((D**3)/(2*(beta**3))))