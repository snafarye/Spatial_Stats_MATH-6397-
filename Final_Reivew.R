getwd()
load('~/Grad_School/MATH 6397/FinalReview.RData')
library(fields)
library(geoR)
library(gstat)

#--------------------------------------------------------------------
#                        Data info


# • x1 and y1 contain location information for the data
# locations to be used for model fitting
# • corresponding spatial variable as well as covariate
# are given in z1 and w1, respectively
# • x0 and y0 give Kriging location
# • z0 and w0 for the spatial variable and covariate
# values for those locations
# 
loc   = cbind(x1,y1)
loc_k = cbind(x0,y0)
quilt.plot(loc, z1, xlim=c(-1,11), ylim=c(-1,3.5))
world(add=T)
points(loc_k, col="black", lwd = 3)  # add in where the kriged points are located

# z0 values at those 5 kriged locations 
quilt.plot(loc_k, z0, xlim=c(-1,11), ylim=c(-1,3.5))
world(add=T)

#--------------------------------------------------------------------
#                        Isotropic modeling with LS
# Part1: expo cov function, OLS

vario0  = variog(coords=loc,
                 data=z1,
                 trend = "2nd"
                 ,max.dist = 3  # can set a diff max distance 
                 ) 
plot(vario0, main = "Emperical Variogram")

alpha  = 6
beta   = 0.25
nugget = 2
# exponential variogram
plot(vario0, main = "EXPO" ) 
lines.variomodel(cov.model = "exp", 
                 cov.pars = c(alpha,beta), 
                 nugget = nugget, 
                 max.dist = 4,
                 lwd = 3)

# exponential covariance parameters 
fit.expo.OLS = variofit(vario0, ini.cov.pars=c(alpha,beta),
                      weights = 'equal',
                      cov.model = "exponential")
fit.expo.OLS



##################
#  Part2: Gaussian cov function, OLS

alpha.g = 6
beta.g  = 0.25
nugget.g = 2

# Gaussian variogram
plot(vario0, main = "Gaussian" ) 
lines.variomodel(cov.model = "gaussian", 
                 cov.pars = c(alpha.g,beta.g), 
                 nugget = nugget.g , 
                 max.dist = 4,
                 lwd = 3)

# Gaussian covariance parameters 
fit.gaussian.OLS = variofit(vario0, ini.cov.pars=c(alpha.g,beta.g),
                        weights = 'equal',
                        cov.model = "gaussian")
fit.gaussian.OLS


#--------------------------------------------------------------------
#--------------------------------------------------------------------
#                        Isotropic modeling with REML
# Part1: expo cov function, OLS

fit.exp.reml = likfit(coords=loc, 
                    data=z1,
                    lik.method = "REML", 
                    trend = "2nd",
                    cov.model = "exponential",
                    ini.cov.pars= c(alpha,beta)) 
fit.exp.reml





fit.gaussian.reml = likfit(coords=loc, 
                      data=z1,
                      lik.method = "REML", 
                      trend = "2nd",
                      cov.model = "gaussian",
                      ini.cov.pars= c(alpha.g,beta.g)) 
fit.gaussian.reml


#--------------------------------------------------------------------
#--------------------------------------------------------------------
#                       plot the variogram fits

par(mfrow=c(1,2))
plot(vario0, main = "EXPO vs Gaussian OLS")

lines(fit.expo.OLS,      lwd = 2, col = "green")
lines(fit.gaussian.OLS,  lwd = 2, col= "red")
legend("bottomright", legend=c("expo.OLS", "gaussian.REML"),
       lwd=c(2,2,2),col = c("green","red"), cex=0.7)

plot(vario0, main = "EXPO vs Gaussian REML")
lines(fit.exp.reml,      lwd = 2, col= "blue")
lines(fit.gaussian.reml, lwd = 2, col= "orange")
legend("bottomright", legend=c("expo.REML", "Gaussian.REML"),
       lwd=c(2,2,2),col = c("blue", "orange"), cex=0.7)


par(mfrow=c(1,1))


print(list(fit.expo.OLS,
             fit.gaussian.OLS,
             fit.exp.reml,
             fit.gaussian.reml)
      )



#--------------------------------------------------------------------
#--------------------------------------------------------------------
#                        Local stationary model

{
data = as.data.frame(cbind(x1,y1,z1,w1))
data1 = data[data$y1 >  1,] 
data2 = data[data$y1 <= 1,] 
loc   = cbind(x1,y1)
quilt.plot(loc, z1, xlim=c(-1,11), ylim=c(-1,3.5))
world(add=T)
points(loc_k, col="black", lwd = 3)  # add in where the kriged points are located
abline(h = 1, col="gray") # ablines at 1 latitude
text(c(0,0),
     c(3,-0.5),
     label= c(nrow(data1), nrow(data2)), 
     col="black") 

text( c(0,0),
      c(3.25,-.25),
      label = c("data1", "data2" ), 
      col = "black")
}


{
vario1=variog(coords=cbind(data1$x1, data1$y1), data=log(data1$z1), trend="2nd",max.dist=10) 
vario2=variog(coords=cbind(data2$x1, data2$y1), data=log(data2$z1), trend="2nd",max.dist=10) 
par(mfrow=c(1,2)) 
plot(vario1,ylim=c(0,0.10),main = "Data1") 
plot(vario2,ylim=c(0,0.10), main= "Data2") 
}














