library("fields") 
library(ggplot2)
library(cowplot)

load("lec10-ex.RData") 


z1
z2
par(mfrow=c(1,2))
plot(z1);plot(z2)
par(mfrow=c(1,1))
quilt.plot(x,y, z1)
quilt.plot(x,y, z2)


############################
#######################
############

#----- using fields package to do the OLS ------------

vg1.1=vgram(cbind(x,y), z1) 
plot(vg1.1) 
# cons, variogram dec as distance is inc


vg1.2=vgram(cbind(x,y), z1, dmax=0.6) 
plot(vg1.2) 


# take the regression then variogram of the residual, 
# does it make a difference,,,,,,no
reg1=lm(z1~1) 
res1= reg1$residuals
#res1=z1-mean(z1) 
vg1.1.r=vgram(cbind(x,y),res1) 
plot(vg1.1.r)   #----- get the same plot as before -------

par(mfrow = c(1, 2))
plot(vg1.1)
plot(vg1.1.r)
par(mfrow = c(1, 1))

#--- take the dmax, get the same results as vg1.2
vg1.2.r=vgram(cbind(x,y), res1, dmax=0.6) 
plot(vg1.2.r) 



####
#----z2 = B0 + B1*X + B2*Y +e    ------------------------
vg2.1=vgram(cbind(x,y), z2) 
plot(vg2.1)  #---- constantly going up, never stabilizing
# mean is spacially variant -------

vg2.2=vgram(cbind(x,y), z2, dmax=0.6) 
plot(vg2.2) 
# posible spcial vareant mean

par(mfrow = c(1, 2))
plot(vg2.1)
plot(vg2.2)
par(mfrow = c(1, 1))




reg2=lm(z2~1+x+y) 
summary(reg2) # see that xx,and,y are sig covariates 

res2=reg2$residuals 
vg2.1.r=vgram(cbind(x,y),res2) 
plot(vg2.1.r)  # smaller sill, spatial range, nicer pattern
#    looks like before with the constant mean
vg2.2.r=vgram(cbind(x,y), res2, dmax=0.6) 
plot(vg2.2.r) 

par(mfrow = c(1, 2))
plot(vg2.1)
plot(vg2.2.r)
par(mfrow = c(1, 1))


###################
###########
##
#----- using gstats package to do the OLS ------------

install.packages(c("gstat","sp"))

library(gstat) 
library(sp) 

data=data.frame(cbind(x,y,z1,z2)) 
names(data)=c("x","y","z1","z2") # rename the cols

coordinates(data)= ~ x+y 

bubble(data, zcol="z1", fill=TRUE, do.sqrt=FALSE) 
bubble(data, zcol="z2", fill=TRUE, do.sqrt=FALSE) 

vg1.3=variogram(z1~1,data=data) 
plot(vg1.3) # automatically decides the best dmax and 
# how many bins to use

vg2.3=variogram(z2~x+y, data=data) 
plot(vg2.3) 

#---- How to estimate parameters  using least squares 

##"Exp" (for exponential) 
##"Sph" (for spherical) 
##"Gau" (for Gaussian) 
##"Mat" (for Matern) 

model1=vgm(psill=1, model="Exp", nugget=0.001, range=0.2) 
# not to estimate initial guesses, watch video--------
plot(vg1.3,model=model1) 

fit1=fit.variogram(vg1.3, model=model1) 
#    fit1= the answers 
plot(vg1.3, model=fit1) 



vg1.4=variogram(z1~1, data=data, cutoff=0.8, width=0.8/20)  
plot(vg1.4)
model2=vgm(psill=0.8, model="Exp", nugget=0.001, range=0.4) 
plot(vg1.4,model=model2) 
fit2=fit.variogram(vg1.4, model=model1)  
plot(vg1.4, model=fit2) 


##https://gsp.humboldt.edu/olm/R/04_01_Variograms.html 


## ---- excercise --------


## estimate a model for z2 = ~ x+y 
# plot(variogram(z2~x+y, data=data)) 
vg5 = variogram(z2~x+y, data=data)  
plot(vg5)

model3 = vgm(psill = 1, 
             model="Exp", 
             nugget=0.01, 
             range=0.3) 
plot(vg5, model=model3)


fit5 = fit.variogram(vg5, model=model3)  
fit5
plot(vg5, model=fit5)

