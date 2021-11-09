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
# expodentail 
plot(vario1, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
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
plot(vario2, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
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
plot(vario3, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
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
plot(vario4, main = "EXP" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.05,0.05), 
                 nugget = 0.005, 
                 max.dist = 1,
                 lwd = 3)



#---plot-----------------------------------------------------

par(mfrow=c(2,2)) 
plot(vario1, main = "Data1")
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.04,0.05), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)

plot(vario2, main= "Data2") 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.10,0.05), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)
plot(vario3, main= "Data3")
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.10,0.10), 
                 nugget = 0.0005, 
                 max.dist = 1,
                 lwd = 3)
plot(vario4, main= "Data4") 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(0.05,0.05), 
                 nugget = 0.005, 
                 max.dist = 1,
                 lwd = 3)
par(mfrow=c(1,1))


###------ variofit() of the subgroups-------------------------------------

# might want to change the initial parameters ......

fit1=variofit(vario1, ini.cov.pars=c(0.04,0.1),
                      weights = 'equal',
                      cov.model = "exponential")
fit1


fit2=variofit(vario2, ini.cov.pars=c(0.1,0.1),
              weights = 'equal',
              cov.model = "exponential")
fit2


fit3=variofit(vario3, ini.cov.pars=c(0.1,0.1),
              weights = 'equal',
              cov.model = "exponential")
fit3  # no spatial dependence? weak spatial dependence 

fit4=variofit(vario2, ini.cov.pars=c(0.05,0.05),
              weights = 'equal',
              cov.model = "exponential")
fit4

print(list(fit1, fit2, fit3, fit4))



#------- kriging the local stay 
