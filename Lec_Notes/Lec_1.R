#Download “data-lecture1.RData” from TEAMS. The data
#set consists of x, y, and two more variables A and B at
#each spatial location. Using R and software fields, do the
#following:
#  1. Draw a scatter plot of A (x-axis)and B (y-axis).
#2. Perform linear regression on B using A. Summarize your
#results of the regression and discuss the relationship
#between the two variables A and B.
#3. Plot a ”map” of B (use ”quilt.plot” command in R
#                      package fields)
#4. Repeat 1 for the residuals from the regression in 2.
#Discuss similarities or differences betwen the two
#spatial maps of B and the residuals.

#A = sample(x = c(1:10), size = 100, replace = TRUE)
#B = sample(x = c(1:20), size = 100, replace = TRUE)
A
B
x
y


Data = data.frame(A,B)

# part 1
#pairs(Data)
plot(A,B, main = "Scatter plot of A and B")

# part 2 
data.lm = lm(B~A, data = Data)
summary(data.lm)

# y = 0.35565Bo +0.82634B1 +ei

# part 3
plot(x,y)
install.packages("fields")
library(fields)
?quilt.plot
# spaial map of B
data.quilt.plot = quilt.plot(x,y, B) # can also look into A



# part 4  map of the residual
plot(A,B)
data.lm$residuals

quilt.plot(x,y, data.lm$residuals)
# want to take a look at the residual 


par(mfrow=c(1,2))
data.quilt.plot = quilt.plot(x,y, B)
quilt.plot(x,y, data.lm$residuals)
par(mfrow=c(1,1))




















