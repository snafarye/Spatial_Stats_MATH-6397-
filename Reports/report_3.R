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

loc = cbind(sample_data$longitude, sample_data$latitude)
y      <- log(sample_data$median_house_value) #value at loc
y_res  <- res0
x1     <- sample_data$total_bedrooms 
x2     <- sample_data$median_income  
x3     <- sample_data$housing_median_age 
x4     <- sample_data$longitude
x5     <- sample_data$latitude

quilt.plot(sample_data$longitude, sample_data$latitude, y)
US(add=T) 


#----- linear regression model from report 1
fit0 = lm(y ~ x1+x2+x3+x4+x5,  data = sample_data)
summary(fit0)
res0 = fit0$residuals
summary(res0)

#plot(res0)
quilt.plot(sample_data$longitude, sample_data$latitude, res0)
US(add=T) 

#------------------- variogram----------------
vario0 = variog(coords=loc,
                     data=y, 
                     trend=trend.spatial(~x1+x2+x3+x4+x5)) 
plot(vario0) 

# ------------------ OLS ------------------------------------------------
zz=rnorm(dim(sample_data)[1], 0,0.1) 
sample_data$longitudee=sample_data$longitude+zz 
locc = cbind(sample_data$longitudee, sample_data$latitude)


fit.expo.OLS=variofit(vario0, ini.cov.pars=c(0.15,0.2),
                      weights = 'equal',
                      cov.model = "exponential")
fit.expo.OLS
#variofit: model parameters estimated by OLS (ordinary least squares):
#  covariance model is: exponential
#parameter estimates:
#  tausq sigmasq     phi 
#0.0479  0.0542  0.0186 
#Practical Range with cor=0.05 for asymptotic range: 0.05568734
#variofit: minimised sum of squares = 0.0017

D <- rdist.earth(locc,miles = T) # distance

alpha = 0.0542
beta  = 0.0186
delta = 0.0479

M <- cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) # design matrix
S <- alpha*exp(-D/beta)                       # covareiance matrix 
diag(S) = diag(S) + delta                     # nugget
Z = matrix(res0, ncol = 1)

B_estimate = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
B_estimate
print(round(B_estimate,5))

-1.045835e+02
x1  2.060834e-04
x2  1.529699e-01
x3  2.984820e-03
x4 -1.168275e+00
x5 -7.078880e-01

# with res0
3.412858e-02
x1  8.190126e-07
x2 -2.121486e-04
x3  5.032516e-05
x4  6.013972e-04
x5  1.008981e-03

##---------- Kriging  --------------------------------------------------------- 
res <- (lm(log(median_house_value) ~ total_bedrooms+ 
             median_income+housing_median_age + 
             longitude+latitude,  data = sample_data))$residuals

# 100 kirging data
yk      <- log(sample_kriging_data$median_house_value) #value at loc
k1     <- sample_kriging_data$total_bedrooms
k2     <- sample_kriging_data$median_income  
k3     <- sample_kriging_data$housing_median_age 
k4     <- sample_kriging_data$longitude
k5     <- sample_kriging_data$latitude

# locations for selected data and kriging data
locc = cbind(sample_data$longitude, sample_data$latitude)
locc_k = cbind(sample_kriging_data$longitude, sample_kriging_data$latitude)

# distance
D <- rdist.earth(locc,miles = T) # distance b/t the selected data
d = rdist.earth(locc,  locc_k, miles = T)         # distance b/t selected data and kriging data

# modify code below  
M = cbind(rep(1, dim(D)[1]), x1,x2,x3,x4,x5) 
m = t(cbind(rep(1,100), k1,k2,k3,k4,k5)) 

Z = matrix(res, ncol = 1)       # data vector

# -------------------------(my results)
######## exponential OLS
alpha = 0.0542
beta  = 0.0186
delta = 0.0479

S = alpha * exp(-D/beta)         # covariance matrix 
diag(S) = diag(S)+delta          # add nugget 
k = alpha*exp(-d/beta)
lambda = (solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
dim(lambda)

krig.exp.ols = t(lambda) %*% Z 
dim(krig.exp.ols)


plot(yk,krig.exp.ols)

mse1=mean((krig.exp.ols-yk)^2)
mse1 #in log form
mae1=mean(abs(krig.exp.ols-yk))
mae1
lines(yk,yk) #y=x line not perfect predictions but good




# with yk log
> mse1=mean((krig.exp.ols-yk)^2)
> mse1 #in log form
[1] 0.04423849
> mae1=mean(abs(krig.exp.ols-yk))
> mae1
[1] 0.1524849


# with res of the kriging
> mse1=mean((krig.exp.ols-yk)^2)
> mse1 #in log form
[1] 152.5091
> mae1=mean(abs(krig.exp.ols-yk))
> mae1
[1] 12.34164



##-----------  local statinarity ----------------------------

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
data3$total_bedrooms  + 
  data3$median_income +
  data3$housing_median_age +
  data3$longitude+ 
  data3$latitude

vario1=variog(coords=cbind(data1$longitude, data1$latitude), 
              data=log(data1$median_house_value), 
              trend=trend.spatial(~data1$total_bedrooms  + 
                                    data1$median_income +
                                    data1$housing_median_age +
                                    data1$longitude+ 
                                    data1$latitude
                                  )) 
plot(vario1)

vario2=variog(coords=cbind(data2$longitude, data2$latitude), 
              data=log(data2$median_house_value), 
              trend=trend.spatial(~data2$total_bedrooms  + 
                                    data2$median_income +
                                    data2$housing_median_age +
                                    data2$longitude+ 
                                    data2$latitude
                                  )) 
plot(vario2)

vario3=variog(coords=cbind(data3$longitude, data3$latitude), 
              data=log(data3$median_house_value), 
              trend=trend.spatial(~data3$total_bedrooms  + 
                                    data3$median_income +
                                    data3$housing_median_age +
                                    data3$longitude+ 
                                    data3$latitude
                                  ))
plot(vario3)

vario4=variog(coords=cbind(data4$longitude, data4$latitude), 
              data=log(data4$median_house_value), 
              trend=trend.spatial(~data4$total_bedrooms  + 
                                    data4$median_income +
                                    data4$housing_median_age +
                                    data4$longitude+ 
                                    data4$latitude
                                  )) 
plot(vario4)


par(mfrow=c(1,4),mai=c(0.5,0.5,0.5,0.5)) 
plot(vario1,ylim=c(0,1.1),main = "Data1") 
plot(vario2,ylim=c(0,1.1), main= "Data2") 
plot(vario3,ylim=c(0,1.1), main= "Data3") 
plot(vario4,ylim=c(0,1.1), main= "Data4") 
par(mfrow=c(1,1))













