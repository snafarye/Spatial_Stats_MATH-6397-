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
# not much trend, or change in variance 
boxplot(wind1~loc[,1])
#interquartile range, the box, does not seem to change much over time 
boxplot(wind1~loc[,2])
#spatial and temporal mean
sp.mean   = colSums(wind1)/12    # avg space over 12 stations 
time.mean = rowSums(wind1)/365 # avg of time at each 12 station
# the average over the 12 locations 

#no temporal trend or inc/dec pattern
plot(1:365, sp.mean, type="l")  # plot for spatial , using the mean of space

library(viridis)
cols=inferno(100)[100:1]       # plot of the time , using the mean of time
quilt.plot(loc, time.mean, xlim=c(-10.5,-5), ylim=c(51, 56),col=cols)
world(add=T)



# z=matrix(wind1, ncol=1, byrow=FALSE)
# hist(z)
# qqnorm(z);qqline(z, col = "red")
# 
# hist(log(z))
# qqnorm(log(z)); qqline(log(z), col = "red")




##    Seperalble model? 
#-----------------------------------------------------------------------------
# new space time locations , long, lat repeated 365 times 

time=matrix(sort(rep(1:365, 12)),ncol=1)

loc1=NULL
for(i in 1:365){
  loc1=rbind(loc1, loc)
}
dim(loc1) # 4380 = 12*365

# spatrial loc distances 
slag=rdist.earth(loc1, miles=FALSE) ## spatial distance (large)
tlag=rdist(time) ## temporal distance (large)
dim(slag);dim(tlag)



#----need to create a variogram, we have non-constant mean so fit the model + residuals 
## linear regression
# z=matrix(wind1, ncol=1, byrow=FALSE)
# z.mean = matrix(wind.loc$MeanWind, ncol = 1, byrow=FALSE)
# 
# mean1=NULL
# for(i in 1:365){
#   mean1=rbind(mean1, z.mean)
# }
# dim(mean1) # 4380 = 12*365
# 
# 
# 
# fit=lm(mean1~1)
# summary(fit)
## squared difference residuals
# dd=rdist(fit$residuals)
# dd=dd^2; dim(dd)
# dd[1:5,1:5]


# need for the covariance structure, z1 is mean zero 
z=matrix(wind1, ncol=1, byrow=FALSE)
z1=z-mean(z)
dd <- rdist(z1)
dim(dd)





## temporal bin
range(tlag)  # 0 , 364
tbin = 0:12
#seq(0,200,,14)

# %%%%  =========================   time 45min in lec nov 10th

## spatial bin
range(slag) # 0, 427.8
sbin=seq(0,200,,16)


dim(tlag); length(tbin);dim(slag);length(sbin)
dim(dd)



st.vario = array(NA, dim=c(length(tbin)-1, length(sbin)-1))
dim(st.vario)

for(i in 2:length(tbin)){
  for(j in 2:length(sbin)){
    print(i)
    print(j)
    
    temp = dd[slag>=sbin[j-1] & slag<sbin[j] & tlag>=tbin[i-1] & tlag<tbin[i]]
    
    a=mean(temp, na.rm = T)  # space/ time variogram value
    print(a)
    
    st.vario[i-1,j-1] = a
  }
}
# time bin, spatial bin, space/time variogram value 
st.vario[1:5,1:5]
# space time values for computation of spatial lag and time lag 




# image plot of values above 
library(RColorBrewer)

col=brewer.pal(9, "Blues")
image.plot(tbin[1:length(tbin)], 
           sbin[1:length(sbin)],
           st.vario,
           xlab="temporal lag (month)",
           ylab="spatial lag (mile)",
           col=col)


plot(tbin[1:length(tbin)-1], st.vario[,1]) # 
plot(sbin[1:length(sbin)-1], st.vario[1,]) #
# sign of sig spatial dependence ????????



#
plot(sbin[1:length(sbin)-1], seq(0, 1.2,,15), 
     type="n", 
     xlab="spatial lag",
     ylab="variogram")
col1=brewer.pal(12, "Paired")
for(i in 1:12){
  points(sbin[1:length(sbin)-1],
         st.vario[i,],
         col=col1[i],pch=20)
}

# spatial variogram , diff color = diff temporal lag 




#
plot(tbin[1:length(tbin)-1], 
     seq(0, 1.2,, 12), 
     type="n", 
     xlab="temporal lag",
     ylab="variogram")
col1=brewer.pal(12, "Paired")
for(i in 1:12){
  points(tbin[1:length(tbin)-1],
         st.vario[,i],
         col=col1[i],
         pch=20)
}
# temporal variogram , diff color = diff spatial lag 
# temporal dependence is weak (compared to spatial dependence )




