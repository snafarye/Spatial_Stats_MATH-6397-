##  lec 11
load('lec10-ex.RData')

install.packages('geoR')
library(geoR)





#---  z1 = B0 + e       (unknonw constant mean and cov) 



# sigma = alpha
# phi = spacial range
# tausq = nugget 
# beta = mean parameter 
# -------------------------------------------------------------------------------------
# -------------   want a exp w/nuggest
fit1 = likfit(coords=cbind(x,y), data  = z1, 
              fix.nugget = TRUE, nugget = 0,  # want the nugget to not be estimated, but fixed at 0 (ie. no nugget) 
              fix.kappa = TRUE, kappa = 0.5, #-- exponential example
              ini.cov.pars = c(1,0.3))       # initial values for alpha = 1, and range = 0.3
fit1

# cov.model at default matern

#     beta      sigmasq      phi 
#    "0.8165"   "0.7029"    "0.1872" 
#Practical Range with cor=0.05 for asymptotic range: 0.5609485
#likfit: maximised log-likelihood = -78.19





# -------------------------------------------------------------------------------------
# -----------    want a exp w nugget
fit2 = likfit(coords=cbind(x,y), data  = z1, 
              fix.nugget = FALSE,              # have it estimate a nugget 
              fix.kappa = TRUE, kappa = 0.5,  #-- exponential example
              ini.cov.pars = c(1,0.3))         # initial values for alpha = 1, and range = 0.3
fit2  
# likfit: maximised log-likelihood = -75.16
#     beta    tausq  sigmasq      phi 
#   "0.8757" "0.0305" "0.7207" "0.2716" 
# Practical Range with cor=0.05 for asymptotic range: 0.8135489


# get a higher value  maximize log-likelihood than fit1, (and that's why you want there for its better)
# dif, cov matrix changes, the pdf is the same tho
# but cant just say the second one is better, can use AIC for more accuracy, 
# need to have same distribution, Mulitval distribution, 
# but comparison not valid fit1 and fit2 have diff parameters
# we should look at AIC instead of just the log-likelihood



# -------------------------------------------------------------------------------------
# -------       want estimate smoothness
fit3 = likfit(coords=cbind(x,y), data  = z1, 
              fix.nugget = FALSE,
              fix.kappa = FALSE,       # make it estimate kappa
              ini.cov.pars = c(1,0.3))  # initial values for alpha = 1, and range = 0.3
fit3
#     beta    tausq        sigmasq      phi       kappa 
# "0.8757"    "0.0305"    "0.7207"    "0.2716"    "0.5000" 
# Practical Range with cor=0.05 for asymptotic range: 0.8135489

# likfit: maximised log-likelihood = -75.16


# no diff from fit2, since default is kappa = 0.5 --> exponential 
# -------------------------------------------------------------------------------------

fit4 = likfit(coords=cbind(x,y), data  = z1, 
              fix.nugget = FALSE,        # estimate a nugget
              fix.kappa = FALSE,
              ini.cov.pars = c(1,0.3),   # initial values for alpha = 1, and range = 0.3
              cov.model = 'exponential')        # cov. fun. model = : "matern", "exponential", "gaussian", "wave"
fit4




# -------------------------------------------------------------------------------------

## ---------   z2  = B0 + B1*X + B2*Y +e

fit5 = likfit(coords=cbind(x,y), data  = z2, 
              fix.nugget = FALSE,
              fix.kappa = TRUE, trend = '1st',   # "1st"the mean is assumed to be a first order polynomial on the coordinates:
                                                    # mu(x) = beta0 + beta1*x1 + beta2*x2.
              ini.cov.pars = c(1,0.3))
fit5
#   beta0    beta1       beta2      tausq       sigmasq      phi 
# "0.2194"   "3.7404"    "2.0905"   "0.0304"    "0.6317"     "0.2331" 
# Practical Range with cor=0.05 for asymptotic range: 0.6982262

# likfit: maximised log-likelihood = -74.78
