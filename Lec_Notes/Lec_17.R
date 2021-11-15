##install.packages("spacetime")
##library(spacetime)

library(gstat)
library(sp)
library(fields)


data("wind") # in the gstat package 
head(wind)
# change long, lat to numerical values 
wind.loc$y <- as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
wind.loc$x <- as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))

plot(wind.loc$x, wind.loc$y,xlim=c(-10.5,-5), ylim=c(51, 56),pch=20, col=2)
world(add=T)
text(wind.loc$x, wind.loc$y, label=wind.loc$Code,col="purple")

#combine time into 1 column
wind$time <- ISOdate(wind$year + 1900, wind$month, wind$day)
wind$jday <- as.numeric(format(wind$time, "%j"))

wind1=wind[1:365, c(4:15)]
wind1=t(wind1)
loc=cbind(wind.loc$x, wind.loc$y) ## longitude and latitude

plot(1:365, seq(0, max(wind1), ,365), type="n")
# time series , over time, diff color = diff stations 
for(i in 1:12){
  lines(wind1[i,],col=i)
}
# is there a pattern b/c long, lat, does the wind speed change with long, lat, each box for each loc, y = speed
# not much treand, or change in vaeriance 
boxplot(wind1~loc[,1])

#interquartile range, the box, does not seem to change much over time 
boxplot(wind1~loc[,2])

#spatial and temporal mean
sp.mean=colSums(wind1)/12  #avg space over 12 stations 
time.mean=rowSums(wind1)/365 # avg of time at each 12 station

# the average over the 12 locations 
#no temporal trend or inc/dec pattern

plot(1:365, sp.mean,type="l")  #time plot

library(viridis)

cols=inferno(100)[100:1]       # space plot
quilt.plot(loc, time.mean,xlim=c(-10.5,-5), ylim=c(51, 56),col=cols)
world(add=T)


z=matrix(wind1, ncol=1, byrow=FALSE)
hist(z)
qqnorm(z);qqline(z, col = "red")
hist(log(z))
qqnorm(log(z)); qqline(log(z), col = "red")







#
# new space time locations , long, lat repeated 365 times 

time=matrix(sort(rep(1:365, 12)),ncol=1)

loc1=NULL
for(i in 1:365){
  loc1=rbind(loc1, loc)
}
dim(loc1) # 4380 = 12*365

# spatrial loc distances 
D.s=rdist.earth(loc1, miles=FALSE) ## spatial distance (large)
D.t=rdist(time) ## temporal distance (large)

D.s1=rdist.earth(loc, miles=FALSE) ## spatial distance (small)
D.t1=rdist(1:365)## temporal distance (small)

z1=z-mean(z)


# calc log-likelihood under seperable 
loglik.sep=function(par){
  
  print(par)
  # take the exp to have positive values in the end
  alpha=exp(par[1])
  beta1=exp(par[2])
  beta2=exp(par[3])
  
  S1=exp(-D.s1/beta1)
  S2=exp(-D.t1/beta2)
  S=alpha*kronecker(S1, S2) # s, t under seperable 
  
  part1=determinant(S, logarithm=T)$modulus  # give the log of the determinate 
  part2=t(z1) %*% solve(S) %*% z1
  
  #K=chol(S)
  #part1=sum(log(diag(K)))*2
  
  #temp=forwardsolve(t(K),z1 )
  #temp=backsolve(K, temp)
  #part2=t(z1) %*% temp
  
  temp=part1+part2
  
  print(-temp)
  
  return(temp)
  
}


# 1st approch 
loglik.sym=function(par){
  
  print(par)
  
  alpha=exp(par[1])
  beta1=exp(par[2])
  beta2=exp(par[3])
  
  D=sqrt((D.s/beta1)^2+(D.t/beta2)^2)
    S=alpha*exp(-D)
  
  #K=chol(S)
  #part1=sum(log(diag(K)))*2
    #temp=forwardsolve(t(K),z1 )
  #temp=backsolve(K, temp)
  #part2=t(z1) %*% temp
  
  temp=part1+part2
    print(-temp)
  return(temp)
  
}

ini=c(log(25), log(100), log(100))

fit.sep=optim(ini, loglik.sep)


## results from optim below

> fit.sep
$par
[1] 3.2418530 1.9740379 0.8465352

$value
[1] 16166.71

$counts
function gradient 
112       NA 

$convergence
[1] 0

$message
NULL
##


fit.sym=optim(ini, loglik.sym)

## results from optim below

> fit.sym
$par
[1] 4.141767 6.374203 1.492660

$value
[1] 15229.56

$counts
function gradient 
78       NA 

$convergence
[1] 0

$message
NULL
##

### Kriging

x0=seq(-10.5,-5,,50)
y0=seq(51, 56,,50)

loc0=make.surface.grid(list(x0,y0))
t0=366

## separable model

par=c(3.2418530, 1.9740379, 0.8465352)

alpha=exp(par[1])
beta1=exp(par[2])
beta2=exp(par[3])

S1=exp(-D.s1/beta1)
S2=exp(-D.t1/beta2)

S=alpha*kronecker(S1, S2)

d.s=rdist.earth(loc, loc0, miles=FALSE)
d.t=rdist(1:365, t0)

k=alpha*kronecker(exp(-d.s/beta1), exp(-d.t/beta2))

lambda=solve(S) %*%k

## symmetric model

par=c(4.141767, 6.374203, 1.492660)

alpha=exp(par[1])
beta1=exp(par[2])
beta2=exp(par[3])

D=sqrt((D.s/beta1)^2+(D.t/beta2)^2)

S=alpha*exp(-D)

d.s=rdist.earth(loc1, loc0, miles=FALSE)
d.t=rdist(time, rep(t0, dim(loc0)[1]) )

d=sqrt((d.s/beta1)^2+ (d.t/beta2)^2)
k=alpha*exp(-d)

lambda=solve(S) %*%k
