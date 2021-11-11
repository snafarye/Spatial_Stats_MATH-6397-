library(readxl)
#load('C:/Users/nsara/Documents/Grad_School/MATH 6397/report/data1.RData')
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

# linear regression model from report 1
fit0 = lm(log(median_house_value) ~ total_bedrooms+ 
            median_income+housing_median_age + 
            longitude+latitude,  data = sample_data)
summary(fit0)
res0 = fit0$residuals

##-----------  local statinarity ----------------------------

# quick renaming 
data0  = sample_data
data.v = sample_kriging_data

range(sample_data$longitude)

data1=data0[data0$longitude > -122.204,] 
data2=data0[data0$longitude <= -122.204 & data0$longitude >  -122.397,]
data3=data0[data0$longitude <=  -122.397 & data0$longitude> -122.59,]

par(mfrow=c(1,1),mai=c(0.5,0.5,0.5)) 
quilt.plot(data0$longitude, data0$latitude, res0) 
US(add=T) 
abline(v= c( -122.204, -122.397),col="gray") # ablines at -99 and -97 longatude


text(-122.5,  37.6, label=nrow(data3), col="black") 
text(-122.3,   37.6, label=nrow(data2), col="black") 
text(-122.1,  37.6, label=nrow(data1), col="black") 
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

#---------------------OLS

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

###########beta estimates ##########################

##data 1

set.seed(222)
zz=rnorm(dim(data1)[1], 0, 0.001) 
data1$longitude=data1$longitude+zz 
locc = cbind(data1$longitude, data1$latitude)

D <- rdist(locc) # distance

alpha = 0.0521 #sill
beta  = 0.1048 #range
delta = 0.0341 #nugget?

M <- cbind(rep(1, dim(D)[1]), data1$total_bedrooms,data1$median_income,data1$housing_median_age,data1$longitude,data1$latitude) # design matrix
S <- alpha*exp(-D/beta)                       # covareiance matrix 
diag(S) = diag(S) + delta                     # nugget
Z = matrix(res_krig_data1, ncol = 1)

B_estimate1 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
print(round(B_estimate1,5))

##data 2

set.seed(222)
zz=rnorm(dim(data2)[1], 0, 0.001) 
data2$longitude=data2$longitude+zz 
locc = cbind(data2$longitude, data2$latitude)

D <- rdist(locc) # distance

alpha = 0.0738 #sill
beta  = 0.0172 #range
delta = 0.0349 #nugget?

M <- cbind(rep(1, dim(D)[1]), data2$total_bedrooms,data2$median_income,data2$housing_median_age,data2$longitude,data2$latitude) # design matrix
S <- alpha*exp(-D/beta)                       # covareiance matrix 
diag(S) = diag(S) + delta                     # nugget
Z = matrix(res_krig_data2, ncol = 1)

B_estimate1 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
print(round(B_estimate1,5))

##data 3

set.seed(222)
zz=rnorm(dim(data3)[1], 0, 0.001) 
data3$longitude=data3$longitude+zz 
locc = cbind(data3$longitude, data3$latitude)

D <- rdist(locc) # distance

alpha = 0.0255 #sill
beta  = 0.0145 #range
delta = 0.0555 #nugget?

M <- cbind(rep(1, dim(D)[1]), data3$total_bedrooms,data3$median_income,data3$housing_median_age,data3$longitude,data3$latitude) # design matrix
S <- alpha*exp(-D/beta)                       # covareiance matrix 
diag(S) = diag(S) + delta                     # nugget
Z = matrix(res_krig_data3, ncol = 1)

B_estimate1 = solve(t(M) %*% (solve(S) %*% M) ) %*%  t(M)%*%solve(S) %*% Z
print(round(B_estimate1,5))

#------------------------------ WLS


fit.WLS1=variofit(vario1, ini.cov.pars=c(0.08,0.02),
              cov.model = "exponential")
fit.WLS1

fit.WLS2=variofit(vario2, ini.cov.pars=c(0.1,0.05),
              cov.model = "exponential")
fit.WLS2

fit.WLS3=variofit(vario3, ini.cov.pars=c(0.08,0.05),
              cov.model = "exponential")
fit.WLS3 
