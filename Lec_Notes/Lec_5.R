load("C:/Users/nsara/Documents/Grad_School/MATH 6397/lec5.RData")

install.packages(c("geoR", "fields", "gstat"))

library(geoR)
library(fields)
library(gstat)


#####
####    in class practice 
###### some packages with their functions to plot variograms
# geoR (variog, variog4) 
# fields (vgram)
# gstat (variogram, variogram with "alpha")
# directional vario using ( variog4, variogram with "alpha")


loc     # x, y coordinates 
z3      # spatial variable
   # what they look like
quilt.plot(loc,z3)

#----------  fields package -------------------------
vg1 = vgram(loc, z3)
names(vg1)
# d     --> spatial distance lag
# vgram --> variogram value
# call  --> 
# type  --> 

# useless graph, takes very long time to load, black colud
# plot(vg1$d, vg1$vgram)

# simple variogram curve 
?vgram
plot(vg1)
plot(vgram(loc, z3, 
           id = NULL, 
           d = NULL, 
           lon.lat = FALSE,  # TRUE --> coordinate are in longitude and latitude (not x, y)
           dmax = NULL,      # max distance
           N = NULL,         # how many bins you have
           breaks = NULL,    # where to break the bins are certain points
           type=c("variogram", "covariogram", "correlogram"))
     )
# how it changes when setting the max distance
par(mfrow=c(1,2))
plot(vg1)
plot(vgram(loc, z3, 
           dmax = 0.8,      # max distance
))
par(mfrow=c(1,1))


# boxplot 
boxplotVGram(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), plot=TRUE, plot.args, ...)




#----------  gstat package -------------------------
?variogram



