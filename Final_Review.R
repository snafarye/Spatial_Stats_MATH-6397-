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
LOC = cbind(c(x1,x0),c(y1,y0))
D = rdist(LOC)


loc   = cbind(x1,y1)
loc_k = cbind(x0,y0)
quilt.plot(loc, z1, xlim=c(-1,11), ylim=c(-1,3.5))
world(add=T)
points(loc_k, col="black", lwd = 3)  # add in where the kriged points are located

# z0 values at those 5 kriged locations 
quilt.plot(loc_k, z0, xlim=c(-1,11), ylim=c(-1,3.5))
world(add=T)

#--------------------------------------------------------------------
#                        Isotropic modeling with OLS

ls()
w2 = w1^2
# since we are working with spatial regression structure 
# z = B0 + B1*W + B2*W^2 + e

fit = lm(z1~w1+w2)
res = fit$residuals

vario0 = variog(coords=cbind(x1,y1), 
                data = res, 
                max.dist = 5)
plot(vario0, main = "Emperical Variogram")

# Part1: expo cov function, OLS
fit1.ols = variofit(vario0, max.dist = 5, weights = "equal", cov.model = "exponential")
fit1.ols
#  Part2: Gaussian cov function, OLS
fit2.ols = variofit(vario0, max.dist = 5, weights = "equal", cov.model = "gaussian")
fit2.ols
# gaussian has the better fit since have lower min sum of squares

# emperical variogram of the residual of the linear regression 
plot(vario0, main = "Emperical Variogram")
lines.variomodel(fit1.ols, col = "blue") #using the cov paramerters from fit1.ols.......  
lines.variomodel(fit2.ols, col = "red")
legend("bottomright", legend=c("expo, OLS", "gaussian, OLS"),
       lwd=c(2,2,2),col = c("blue", "red"), cex=0.7)


#--------------------------------------------------------------------
#--------------------------------------------------------------------
#                        Isotropic modeling with REML

# Part1: expo
fit1.reml = likfit(coords = cbind(x1,y1), data = z1, lik.method = "REML", trend = trend.spatial(~w1+w2), ini.cov.pars = c(1,0.1))
fit1.reml

# Part2: gaussian
fit2.reml = likfit(coords = cbind(x1,y1), data = z1, lik.method = "REML", trend = trend.spatial(~w1+w2), ini.cov.pars = c(1,0.1), cov.model = "gaussian")
fit2.reml

plot(vario0, main = "EXPO vs Gaussian REML")
lines(fit1.reml,      lwd = 2, col= "blue")
lines(fit2.reml, lwd = 2, col= "orange")
legend("bottomright", legend=c("expo.REML", "Gaussian.REML"),
       lwd=c(2,2,2),col = c("blue", "orange"), cex=0.7)

#--------------------------------------------------------------------
#--------------------------------------------------------------------
#                       plot the variogram fits

par(mfrow=c(1,2))
plot(vario0, main = "Emperical Variogram")
lines.variomodel(fit1.ols, col = "blue") #using the cov paramerters from fit1.ols.......  
lines.variomodel(fit2.ols, col = "red")
legend("bottomright", legend=c("expo, OLS", "gaussian, OLS"),
       lwd=c(2,2,2),col = c("blue", "red"), cex=0.7)


plot(vario0, main = "EXPO vs Gaussian REML")
lines(fit1.reml,      lwd = 2, col= "blue")
lines(fit2.reml, lwd = 2, col= "orange")
legend("bottomright", legend=c("expo.REML", "Gaussian.REML"),
       lwd=c(2,2,2),col = c("blue", "orange"), cex=0.7)


par(mfrow=c(1,1))

#--------------------------------------------------------------------
fit1.ols
summary(fit)
fit1.reml
#--------------------------------------------------------------------
#                        Local stationary model
loc = cbind(x1,y1)

loc1 = loc[loc[,2] > 1,]
loc2 = loc[loc[,2] <= 1,]

z11 = z1[loc[,2] > 1]
z12 = z1[loc[,2] <= 1]

w11 = w1[loc[,2] > 1]
w12 = w1[loc[,2] <= 1]



w112 = w11^2
fit1.nons = likfit(coords = loc1,
                   data = z11,
                   lik.method = "REML",
                   trend = trend.spatial(~w11+w112),
                   ini.cov.pars = c(1,0.1) )
w122 = w12^2
fit2.nons = likfit(coords = loc2,
                   data = z12,
                   lik.method = "REML", cov.model = "gaussian",
                   trend = trend.spatial(~w12 + w122),
                   ini.cov.pars = c(1,0.1) )


#--------------------------------------------------------------------
#--------------------------------------------------------------------
#       Comparison regarding kriging accuracy 







#--------------------------------------------------------------------
#--------------------------------------------------------------------


#---------------------------------------------------------------------




















