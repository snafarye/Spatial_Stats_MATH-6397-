#install.packages("spatstat")
library(spatstat)


library(readr)
data <- read_csv("spatstat-ex.csv")
head(data)

#convert the data to a point pattern object using the spatstat command ppp
range(data$x)
range(data$y)
mypattern <- ppp(data$x, data$y, c(range(data$x)[1], range(data$x)[2]),
                                   c(range(data$y)[1], range(data$y)[2]) )
plot(mypattern)

##   stats analysis

summary(mypattern)
plot(Kest(mypattern))

plot(envelope(mypattern,Kest))

#This is a spatial map of first order intensity function
plot(density(mypattern))

# add in the other data column 
marks(mypattern) <- data[,c(3,4)]
plot(Smooth(mypattern))












