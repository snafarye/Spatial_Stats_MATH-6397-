set.seed(0920) 
library(fields) 

?Matern

x=runif(100,0,1) 
y=runif(100,0,1) 

#--- calc the distance matrix

D=rdist(cbind(x,y))  # distance already takes into account 
dim(D)
D[1:5,1:5]


##--- exponential with alpha=2, beta=0.2 

S1=2*exp(-D/0.2) 
dim(S1) #covariance matrix
round(S1[1:5,1:5],2) # diagonal are 2
# if covariance model is valid then is should
# positive definite ==> eigenvalues of the covariance matrix be positive 
eigen(S1, only.values = TRUE)

##--- add a nugget delta=0.001 

diag(S1)= diag(S1)+0.001 


##--- matern with nu=1.38(smoothness) , alpha=1.5 , beta=0.1(range)

S2 = 1.5*Matern(D , 
                range = 0.1,
                smoothness = 1.38 )
round(S2[1:5,1:5],2)  #diagonal are 1.5
eigen(S2, only.values = TRUE)


##--- gaussian with alpha=1 , beta=0.6

S3 = 1*exp( - (D/0.6)^2 ) 
round(S3[1:5,1:5],2)  #diagonal are 1
eigen(S3, only.values = TRUE)










##---- Using longitude and Latitude 

lat=runif(100,-40,40) 
lon=runif(100,-80,80) 
# when data is given in log, lat use rdist.eath()
# (miles= T --> miles) if (miles= F --> km)
D=rdist.earth(cbind(lon,lat), miles=F) 
D[1:5,1:5]
range(D)
##- - matern with alpha=1, beta=1000km, nu=0.3 
S4= 1*Matern(D, range = 1000, smoothness=0.3) 
S4[1:5,1:5]





x=runif(100,0,1) 
y=runif(100,0,1) 
D=rdist(cbind(x,y))
S5 = exp(- (D/0.9)^2)
S5[1:5,1:5]
eigen(S5, only.values = T)
#add in a nugget
diag(S5)=diag(S5)+10^(-5)
eigen(S5, only.values = T)
## becareful in adding nugget, 
# now the lowest value you have is ((1.0)*10^(-5))
# make sure model is stable, positive def…

