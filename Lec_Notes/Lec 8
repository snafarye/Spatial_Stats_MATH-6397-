load("~/Grad_School/MATH 6397/lec8.RData")


set.seed(100) 
x=runif(25,0,1) 
y=runif(25,0,1) 
library(fields) 
D=rdist(cbind(x,y)) 
S=exp(-D/0.3) 
diag(S)= diag(S) +0.2 
z=rnorm(25,0,1) 
z1=t(chol(S)) %*% z 
z1=z1+x 
z=c(z1[1:20], rep(NA, 5)) 
save(x,y,z, file="lec8.RData") 


library(fields)
quilt.plot(x,y,z)
points(x[21:25],y[21:25]) # the last 5 values that are missing in z but are in x,y

## no spatial dependence? 

z.krig0=rep(mean(z[1:20]),5) #using the mean of 20 values in z to fill in for that last 5 values in z


## simple Kriging ##  

d0=rdist(cbind(x[1:20], y[1:20]), cbind(x[21:25], y[21:25])) 
D0=rdist(cbind(x[1:20], y[1:20])) 
S=exp(-D0/0.3) 
k=exp(-d0/0.3) 
lambda1= solve(S) %*% k  #-- weight
z.krig1=t(lambda1) %*% z[1:20] + (-1)  


## ordinary Kriging ## 

M=matrix(rep(1,20),ncol=1) 
m=matrix(rep(1,5),nrow=1) 
lambda2=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig2=t(lambda2) %*% z[1:20] 


## universal Kriging  ##

M=cbind(rep(1,20), x[1:20]) 
m=t(cbind(rep(1,5), x[21:25])) 
lambda3=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig3=t(lambda3) %*% z[1:20] 

cbind(z.krig0, z.krig1, z.krig2, z.krig3, z1[21:25]) 


## simple Kriging with nugget ##  

d0=rdist(cbind(x[1:20], y[1:20]), cbind(x[21:25], y[21:25])) 
D0=rdist(cbind(x[1:20], y[1:20])) 
S=exp(-D0/0.3) 
      diag(S)=diag(S)+0.2   #-- adding the 0.2 for nugget 

k=exp(-d0/0.3) 
lambda1= solve(S) %*% k 
z.krig1.1=t(lambda1) %*% z[1:20] + (-1)  



## ordinary Kriging with nugget ## 

M=matrix(rep(1,20),ncol=1) 
m=matrix(rep(1,5),nrow=1) 
lambda2=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig2.1=t(lambda2) %*% z[1:20] 


## universal Kriging  with nugget 

M=cbind(rep(1,20), x[1:20]) 
m=t(cbind(rep(1,5), x[21:25])) 
lambda3=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig3.1=t(lambda3) %*% z[1:20] 


cbind(z.krig0, z.krig1, z.krig1.1,z.krig2, z.krig2.1, z.krig3,z.krig3.1, z1[21:25]) 

#----------------- now wth MSEs ------------------------# 

## simple Kriging with nugget ##  

d0=rdist(cbind(x[1:20], y[1:20]), cbind(x[21:25], y[21:25])) 
D0=rdist(cbind(x[1:20], y[1:20])) 
S=exp(-D0/0.3) 
diag(S)=diag(S)+0.2 #-- add in nugget
k=exp(-d0/0.3) 
lambda1= solve(S) %*% k 
z.krig1.1=t(lambda1) %*% z[1:20] + (-1)  

d00=rdist(cbind(x[21:25], y[21:25])) #--- MSE 
k0=exp(-d00/0.3) 
diag(k0)=diag(k0)+0.2 

g= m-t(M) %*% solve(S) %*% k 
# -- note that k0 is the same for each kriging 
mse1=k0-t(k) %*% solve(S) %*% k + t(g) %*% solve(t(M) %*% solve(S) %*% M) %*% g 
mse1=diag(mse1) 
t(t(mse1))



## ordinary Kriging with nugget ## 

M=matrix(rep(1,20),ncol=1) 
m=matrix(rep(1,5),nrow=1) 
lambda2=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig2.1=t(lambda2) %*% z[1:20] 

g= m-t(M) %*% solve(S) %*% k 
# -- note that k0 is the same for each kriging 
mse2=k0-t(k) %*% solve(S) %*% k + t(g) %*% solve(t(M) %*% solve(S) %*% M) %*% g 
mse2=diag(mse2) 

## universal Kriging  with nugget 

M=cbind(rep(1,20), x[1:20]) 
m=t(cbind(rep(1,5), x[21:25])) 
lambda3=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig3.1=t(lambda3) %*% z[1:20] 

g= m-t(M) %*% solve(S) %*% k 
# -- note that k0 is the same for each kriging 
mse3=k0-t(k) %*% solve(S) %*% k + t(g) %*% solve(t(M) %*% solve(S) %*% M) %*% g 
mse3=diag(mse3) 


cbind(z.krig0, z.krig1, z.krig1.1,z.krig2, z.krig2.1, z.krig3,z.krig3.1, z1[21:25]) 



###------ Kriging over the entire domain --------#### 


x1=seq(0,1,,60) 
y1=seq(0,1,,60) 

loc1=make.surface.grid(list(x1,y1)) 

d0=rdist(cbind(x[1:20], y[1:20]), loc1) 
k=exp(-d0/0.3) 

d00=rdist(loc1) 

k0=exp(-d00/0.3) 
diag(k0)=diag(k0)+0.2 

m=t(cbind(rep(1,3600), loc1[,1])) 

lambda4=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 

z.krig4=t(lambda4) %*% z[1:20] 

g= m-t(M) %*% solve(S) %*% k 

mse4=k0-t(k) %*% solve(S) %*% k + t(g) %*% solve(t(M) %*% solve(S) %*% M) %*% g 
mse4=diag(mse4) 

par(mfrow=c(1,2)) 
##  par(mfrow =c(1,1), mai = c(0.1,0.1,0.1,0.1))
image.plot(as.surface(loc1, z.krig4),zlim=c(-1.45, 2.59)) 
quilt.plot(x[1:20], y[1:20], z[1:20], zlim=c(-1.45, 2.59), add=T,nx=60, ny=60) 


image.plot(as.surface(loc1, mse4)) 

points(x[1:20], y[1:20]) 

#--   try 
#   z = b0 + b1 x + b2 y +e
#   cov structure = 2*exp(-d/0.1) with nugget = 0.1






a = 2
b = 0.1
n = 0.1
d0=rdist(cbind(x[1:20], y[1:20]), cbind(x[21:25], y[21:25])) 
D0=rdist(cbind(x[1:20], y[1:20])) 
S=a*exp(-D0/b) 
diag(S)=diag(S)+n   #-- adding the 0.1 for nugget 

k=a*exp(-d0/b) 

M=cbind(rep(1,20), x[1:20]) 
m=t(cbind(rep(1,5), x[21:25])) 
lambda4=(solve(S) - solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% t(M) %*% solve(S) ) %*% k + solve(S) %*% M %*% solve(t(M) %*% solve(S) %*% M) %*% m 
z.krig4=t(lambda4) %*% z[1:20] 

d00=rdist(cbind(x[21:25], y[21:25])) #--- MSE for the parts we wanting to do kriging prediction
k0=a*exp(-d00/b) 
diag(k0)=diag(k0)+n      #-- adding the 0.1 for nugget

g= m-t(M) %*% solve(S) %*% k 
# -- note that k0 is the same for each kriging 
mse4=k0-t(k) %*% solve(S) %*% k + t(g) %*% solve(t(M) %*% solve(S) %*% M) %*% g 
mse4=diag(mse4) 
t(t(mse4))

cbind(mse1,mse2,mse3,mse4)




