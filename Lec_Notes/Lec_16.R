

load("~/Grad_School/MATH 6397/airbnb-all.RData")


loc=data[,c(6,5)]
data[1:6,1:6]
names(data)
range(data$Year)

# create a column of time coordinates, 
# in tearms of total number of months 
time=(data[,3]-2008)*12 + data[,4]
range(time)

data$time=time
data[1:7,1:7]


bed=data[,2]            # bedrooms_count
lprice = log(data[,1])  # log of average_rate_per_night

plot(bed, lprice)
plot(time, lprice) # non-linear trend, variance inc, b/c have more data points in later time

time1=sort(unique(time))


library("fields")

## lets limit time to 41 to 100

lprice1 = lprice[time >= 41]
bed1    = bed[time    >= 41]
loc1    = loc[time    >= 41,]

time1  =  time[time >= 41]
time2  =  sort(unique(time1))


library(effects)
devAskNewPage(ask=TRUE)

for(i in 1:60){
  temp=lprice1[time1==time2[i]]
  quilt.plot(loc1[time1==time2[i],], 
             temp,
             xlim=c(-103.8, -93.6), 
             ylim=c(  25.8,  35.3),main=time2[i], 
             zlim=c(   2.3,   9.211))
  US(add=T)
}


#library(gstat)
## spatial points
#library(sp)
#data1=data[data$Year>2011,]
#temp=(data1$Year-min(data1$Year))*12+data1$Month
#data1=data1[order(temp),]
#sp=SpatialPoints(cbind(x=data1$longitude, y=data1$latitude))
#time=yearmon(cbind(data1$Year, data1$Month))
#mydata=data.frame(log(data1$average_rate_per_night))
#stfdf=STFDF(sp,time,mydata) ?!
#    vv=variogramST(lprice1~bed1, )

###----------------------------------------------------------------------------
###----------------------------------------------------------------------------

# creating index , years after 2011, so starting at 2012
data1=data[data$Year>2011,]
temp=(data1$Year-min(data1$Year))*12+data1$Month
data1$temp <- temp
data1=data1[order(temp),]
head(data1, n = 7)
range(data1$temp)

## random sample a manageable size of data
index=sample(1:dim(data1)[1], 5000)
data2=data1[index,]

# bin spatial distance and temporal lag, 
# in each bin calc the square diff in the bin and take average  

## time lag calc
tlag=rdist(data2$temp)
dim(tlag)
tlag[1:5,1:5]

## spatial lag calc
slag=rdist.earth(cbind(data2$longitude, data2$latitude))
#spatial distance in miles, off diag= pairwise distance 
slag[1:5,1:5]

#----need to create a variogram, we have non-constant mean so fit the model + residuals 
## linear regression
fit=lm(log(data2$average_rate_per_night)~data2$bedrooms_count)
summary(fit)

## squared difference residuals
dd=rdist(fit$residuals)
dd=dd^2
dd[1:5,1:5]


## temporal bin
tbin = 0:12

# %%%%  =========================   time 45min in lec nov 10th

## spatial bin
sbin=seq(0,100,,15)
st.vario=array(NA, dim=c(length(tbin)-1, length(sbin)-1))
for(i in 2:length(tbin)){
  for(j in 2:length(sbin)){
    print(i)
    print(j)
    
    temp=dd[slag>=sbin[j-1] & slag<sbin[j] & tlag>=tbin[i-1] & tlag<tbin[i]]
    
    a=mean(temp,na.rm=T)
    print(a)
    
    st.vario[i-1,j-1]=a
  }
}
# time bin, spatial bin, space/time variogram value 
st.vario[1:5,1:5]
# space time values for computation of spatial lag and time lag 




# image plot of values above 
library(RColorBrewer)

col=brewer.pal(9, "Blues")
image.plot(tbin[1:12], 
           sbin[1:14],
           st.vario,
           xlab="temporal lag (month)",
           ylab="spatial lag (mile)",
           col=col)

plot(tbin[1:12], st.vario[,1]) # sill = 0.65, range = 1(1month) 
plot(sbin[1:14], st.vario[1,]) # sill = 0.90, range = 30(30miles)


#
plot(sbin[1:14], seq(0, 1.2,,14), 
     type="n", 
     xlab="spatial lag",
     ylab="variogram")
col1=brewer.pal(12, "Paired")
for(i in 1:12){
   points(sbin[1:14],
          st.vario[i,],
          col=col1[i],pch=20)
}



#
plot(sbin[1:14], seq(0, 1.2,,14), 
     type="n", 
     xlab="spatial lag",
     ylab="variogram")

for(i in 1:9){
  points(sbin[1:14],
         st.vario[i,],
         col=col[9-i+1],
         pch=20)
}


#
plot(tbin[1:12], 
     seq(0, 1.2,,12), 
     type="n", 
     xlab="temporal lag",
     ylab="variogram")
col1=brewer.pal(12, "Paired")
for(i in 1:12){
    points(tbin[1:12],
           st.vario[,i],
           col=col1[i],
           pch=20)
}
# temporal dependance is weak (compared to spatial dependance )





