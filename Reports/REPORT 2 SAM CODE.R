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
                data=y, 
                trend=trend.spatial(~x1+x2+x3+x4+x5) ) 
plot(vario1) 



#ESTIMATION METHODS - SPHERICAL OLS
ols <- variofit(vario1,
                weights = "equal",
                fix.kappa=FALSE, 
                cov.model = "spherical",
                ini.cov.pars=c(0.15,0.2)) 


#PARAMETER ESTIMATION FUNCTION FOR OLS

D <- rdist.earth(locc,miles = T) # distance

#COVARIANCE PARAMETER ESTIMATES FROM OLS variofit()
alpha = 0.0576
beta  = 0.0627
delta = 0.0458




#ESTIMATION METHODS - SPHERICAL WLS
wls <- variofit(vario1,
                # weights = "equal",
                fix.kappa=FALSE, 
                cov.model = "spherical",
                ini.cov.pars=c(0.15,0.2)) 


#PARAMETER ESTIMATION FUNCTION FOR WLS

D <- rdist.earth(locc,miles = T) # distance

#COVARIANCE PARAMETER ESTIMATES FROM WLS variofit()
alpha = 0.0643
beta  = 0.0729
delta = 0.0459



#ESTIMATION METHODS - SPHERICAL REML

reml <- likfit(coords= cbind(sample_data$longitudee,sample_data$latitude),
               data = log(sample_data$median_house_value), 
               trend=trend.spatial(~x1+x2+x3+x4+x5), 
               lik.method = "REML",
               cov.model = "spherical",
               ini.cov.pars=c(0.15,0.2)) 
