



##-----------  local statinarity ----------------------------
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
                 nugget = 0.005, 
                 max.dist = 0.6,
                 lwd = 3)

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

#---------------variogram estimation  for data1,2,3-------------------------------
res_krig_data1 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data1))$residuals

vario1=variog(coords=cbind(data1$longitude, data1$latitude), 
              data=res_krig_data1, 
              trend=trend.spatial(~data1$total_bedrooms  + 
                                    data1$median_income +
                                    data1$housing_median_age +
                                    data1$longitude+ 
                                    data1$latitude ), max.dist = 0.6) 


res_krig_data2 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data2))$residuals

vario2=variog(coords=cbind(data2$longitude, data2$latitude), 
              data=res_krig_data2, 
              trend=trend.spatial(~data2$total_bedrooms  + 
                                    data2$median_income +
                                    data2$housing_median_age +
                                    data2$longitude+ 
                                    data2$latitude), max.dist = 0.6)


res_krig_data3 <- (lm(log(median_house_value) ~ total_bedrooms+ 
                        median_income+housing_median_age + 
                        longitude+latitude,  data = data3))$residuals

vario3=variog(coords=cbind(data3$longitude, data3$latitude), 
              data=res_krig_data3, 
              trend=trend.spatial(~data3$total_bedrooms  + 
                                    data3$median_income +
                                    data3$housing_median_age +
                                    data3$longitude+ 
                                    data3$latitude), max.dist = 0.39)

###-----------REML mean and cov paramerters -----------------------------------

fit.exp.reml=function(data, res, par){
  
  zz=rnorm(dim(data)[1], 0,0.1) 
  data$longitudee=data$longitude+zz 
  locc = cbind(data$longitudee, data$latitude)
  
  a=par[1] 
  b=par[2]
 
   fit.ER = likfit(coords=locc, 
                    data=res, 
                    trend=trend.spatial(~data$total_bedrooms+data$median_income+data$housing_median_age+
                                          data$longitude+data$latitude), 
                    lik.method = "REML",
                    cov.model = "exponential",
                    ini.cov.pars=c(a, b)) 
  return(fit.ER)
}

fit.exp.reml(data1, res_krig_data1, c(0.08,0.02)  )
fit.exp.reml(data2, res_krig_data2, c(0.10,0.05)  )
fit.exp.reml(data3, res_krig_data3, c(0.08,0.05)  )


> fit.exp.reml(data1, res_krig_data1, c(0.08,0.02)  )
kappa not used for the exponential correlation function

  likfit: end of numerical maximisation.
likfit: estimated model parameters:
  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
"-5.1308" "-0.0001" "-0.0305" "-0.0001" "-0.0373" " 0.0211" " 0.0449" " 0.0378" " 0.0761" 
Practical Range with cor=0.05 for asymptotic range: 0.2279925

likfit: maximised log-likelihood = -14.45


---------------------------------------------------------------
> fit.exp.reml(data2, res_krig_data2, c(0.10,0.05)  )
kappa not used for the exponential correlation function

  likfit: end of numerical maximisation.
likfit: estimated model parameters:
  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
"51.4883" " 0.0000" "-0.0119" "-0.0026" " 0.4271" " 0.0246" " 0.0691" " 0.0574" " 0.1200" 
Practical Range with cor=0.05 for asymptotic range: 0.3595061

likfit: maximised log-likelihood = -102.9



---------------------------------------------------------------
> fit.exp.reml(data3, res_krig_data3, c(0.08,0.05)  )
kappa not used for the exponential correlation function

  likfit: end of numerical maximisation.
likfit: estimated model parameters:
  beta0      beta1      beta2      beta3      beta4      beta5      tausq    sigmasq        phi 
"-78.4588" "  0.0000" " -0.0051" " -0.0022" " -0.6153" "  0.0820" "  0.0621" "  0.0582" "  0.2797" 
Practical Range with cor=0.05 for asymptotic range: 0.8379358

likfit: maximised log-likelihood = -38.58


#######################################################################
####          sperical
###################################################################


fit.sph.reml=function(data, res, par){
  
  zz=rnorm(dim(data)[1], 0,0.1) 
  data$longitudee=data$longitude+zz 
  locc = cbind(data$longitudee, data$latitude)
  
  a=par[1] 
  b=par[2]
  
  fit.ER = likfit(coords=locc, 
                  data=res, 
                  trend=trend.spatial(~data$total_bedrooms+data$median_income+data$housing_median_age+
                                        data$longitude+data$latitude), 
                  lik.method = "REML",
                  cov.model = "spherical",
                  ini.cov.pars=c(a, b)) 
  return(fit.ER)
}

fit.sph.reml(data1, res_krig_data1, c(0.08,0.02)  )
fit.sph.reml(data2, res_krig_data2, c(0.10,0.10)  )
fit.sph.reml(data3, res_krig_data3, c(0.08,0.10)  )
---------------------------------------------------------------
  
  
  
  
> fit.sph.reml(data1, res_krig_data1, c(0.08,0.02)  )
kappa not used for the spherical correlation function

  likfit: end of numerical maximisation.
likfit: estimated model parameters:
  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
"18.7694" " 0.0000" "-0.0275" "-0.0001" " 0.1766" " 0.0787" " 0.0440" " 0.0451" " 0.2114" 
Practical Range with cor=0.05 for asymptotic range: 0.2113997

likfit: maximised log-likelihood = 1.866


---------------------------------------------------------------
  

> fit.sph.reml(data2, res_krig_data2, c(0.10,0.10)  )
kappa not used for the spherical correlation function

  likfit: end of numerical maximisation.
likfit: estimated model parameters:
  beta0     beta1     beta2     beta3     beta4     beta5     tausq   sigmasq       phi 
" 2.8185" " 0.0000" "-0.0108" "-0.0030" " 0.0131" "-0.0282" " 0.0740" " 0.0327" " 0.1329" 
Practical Range with cor=0.05 for asymptotic range: 0.1329406

likfit: maximised log-likelihood = -111.6



---------------------------------------------------------------
  
> fit.sph.reml(data3, res_krig_data3, c(0.08,0.10)  )

  
  likfit: end of numerical maximisation.
likfit: estimated model parameters:
  beta0      beta1      beta2      beta3      beta4      beta5      tausq    sigmasq        phi 
"-27.2342" "  0.0000" " -0.0036" " -0.0023" " -0.2032" "  0.0628" "  0.0638" "  0.0354" "  0.3607" 
Practical Range with cor=0.05 for asymptotic range: 0.3607032

likfit: maximised log-likelihood = -38.75


