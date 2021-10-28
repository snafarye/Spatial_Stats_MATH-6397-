library(readxl)
load('C:/Users/nsara/Documents/Grad_School/MATH 6397/report/data1.RData')
dim(data1)

#---------------------- Data Cleaning ------------------------------------------

general_data = data1 #Create a dataframe with the loaded dataset
nrow(general_data) #Check the number of rows of the dataset
sum(is.na(general_data)) # have about 20 rows with at least one missing value
colSums(is.na(general_data))  # check to see what columns have the missing values 

new_data <- na.omit(general_data) #create a new dataset without missing values
nrow(new_data) #check number of rows of new dataset without missing values
sum(is.na(new_data)) #check to see if all missing values have been removed

set.seed(111)
sample_data = new_data[sample(1:dim(new_data)[1],1500), ] #create a sample dataset of 1500 rows
nrow(sample_data)#check to see if sample dataset has 1500 rows

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

plot(sample_data$population,           log(sample_data$median_house_value)) # some what of a pos linear trend, but not strong
plot(sample_data$households,           log(sample_data$median_house_value)) # some what of a pos linear trend, but not strong
plot(sample_data$total_rooms,          log(sample_data$median_house_value)) # some what of a pos linear trend, but not strong
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(sample_data$longitude,      log(sample_data$median_house_value)) # not really see a sig. relationship b/t house_median_value
plot(sample_data$latitude,       log(sample_data$median_house_value)) # and long./ lat 
par(mfrow=c(1,1))


# Use the same linear regression model you used in team report 1 (with the same covariates)
# parameters: total_bedrooms,   median_income, 

# linear regression model from report 1
fit0 = lm(log(median_house_value) ~
           total_bedrooms+ 
           median_income, 
          data = sample_data)
summary(fit0)
res0 = fit$residuals
plot(res0)
quilt.plot(sample_data$longitude, sample_data$latitude, res0)


# -----------  part 1-----------------------------------------------------------
# Fill in the following table with estimates using the estimation method 
# mentioned below. (p: the number of covariates.)
