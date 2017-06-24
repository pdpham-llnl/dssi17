###########################
#author: Ana Paula Sales  #
#last edited: 07/06/2016  #
###########################


############################
# v2s
# map from innovation variance to marginal AR variance
############################
v2s = function(v,phi){ return(v /(1-phi^2))} 

############################
# s2v
# map from y variance to innovation variance 
############################
s2v = function(s,phi){ return(s *(1-phi^2))} 

############################
# build.AR.Phi
# Given n (number of time points) and phi
# Returns n-by-n matrix \Phi s a function of \phi as in lecture notes 
############################
build.AR.Phi = function(phi, n){
  sigma = matrix(NA, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in i:n){
      sigma[i,j] = phi^{abs(i-j)}
      sigma[j,i] = phi^{abs(i-j)}
    }
  }
  return(sigma)
}

############################
# AR.sample.mvn
# Given n (nsamples),T (ntimepoints), phi, v
# Simulates n samples from joint distribution x_{1:T} ~ N(0, \Sigma), 
# where \Sigma = s*\Phi, \Phi is a function of \phi as in lecture notes
# Returns n-by-T matrix X
############################
AR.sample.mvn = function(nsamples = 3, #number of samples to draw
                         ntimepoints = 100, #length of time series
                         phi = 0.1, #AR(1) main parameter
                         v =1 # innovation variance
){
  s = v2s(v,phi)
  jointMean = rep(0,ntimepoints)
  jointSigma = s * build.AR.Phi(phi, ntimepoints)
  x = rmvnorm(nsamples, mean = jointMean, sigma = jointSigma)
  return(x)
}

############################
# AR.sample.compositional
# Given n (nsamples),T (ntimepoints), phi, v, x0
# Simulates n samples from the compositional form of joint distribution 
# p(x_{1:T}) = p(x_T|x_{T-1})p(x_{T-1}|x_{T-2})...p(x_2|x_1)p(x_1)  
# Returns n-by-T matrix X
############################
AR.sample.compositional = function(nsamples = 3, #number of samples to draw
                                   ntimepoints = 100, #length of time series
                                   phi = 0.1, #AR(1) main parameter
                                   v = 1, #innovation variance
                                   x0 = 0, #initial value of the series
                                   verbose = FALSE #mostly for debugging purposes
){
  #draw x1 ~ N(0,s)
  #for t = 2:n, sucessively generate (x_t|x_{t-1}): sample epsilon_t ~ N(0,v) and add it to x_{t-1}
  s = v2s(v,phi)
  s.sr = sqrt(s)
  v.sr = sqrt(v)
  X = NULL
  for(ns in 1:nsamples){
    x = c(rnorm(1,0, s.sr) + phi*x0, rnorm(ntimepoints - 1, 0 , v.sr)) #sample independent innovations/noise
    for(i in 2:ntimepoints)
      x[i] = phi*x[i-1] + x[i] #multiply x_{t-1} by phi and add it to epsilon_{t}
    X = rbind(X,x)
  }
  if(verbose)
    print(paste("nsamples:", nsamples, "ntp:", ntimepoints, "phi:", phi, "v:",v, "s:",s, "x0:",x0))
  return(X)
}

############################
# AR.ref.posterior.inf
# Given x_{1:n} ~ AR(1|\phi,v)
# Computes reference posterior parameters as per lecture notes 
# Returns a names list with posterior parameters
############################
AR.ref.posterior.inf = function(x #your AR series observations vector
                                ){
  x = x - mean(x) 
  n = length(x)
  x.tm1 = x[1:(n-1)] #x_{t-1}
  x.t = x[2:n]       #x_t
  B = x.tm1 %*% x.tm1
  b = ( x.tm1 %*% x.t ) / B
  res.vec = x.t - b*x.tm1 #residual
  Qb = res.vec %*% res.vec
  nu = n - 2
  tau = Qb / nu
  return(list(B = B, b = b, res.vec = res.vec, Qb = Qb, nu = nu, tau = tau))
}


############################
# AR.ref.posterior.and.predictive
# 
############################
AR.ref.posterior.and.predictive = function(x, # observations from AR(1) (all we need is the last item, so you can pass a vector or just the last item)
                                           B, #part of posterior variance of phi|x,v
                                           b, #posterior mean of phi |x, v
                                           nu, #number degrees of freedom (for AR(1) it's n-2)
                                           Qb, #residual sum of squares based on fitted parameter value b
                                           nsamples = 10, #number posterior samples
                                           forecast.length = 100 #how many steps into the feature to make predictions for AR process
){
  v.samp = Qb/rchisq(nsamples, df = nu)
  phi.samp = rnorm(nsamples, 0,1) * sqrt(v.samp/B) + b
  x.samp = mapply(AR.sample.compositional, v = v.samp, phi = phi.samp, x0 = x[length(x)], nsamples =1, ntimepoints = forecast.length)
  return(list(vsamp = v.samp, phisamp = phi.samp, xsamp = x.samp))
}



############################
# AR.plot.samples
# Quick and dirty plotting of AR(1) time series
# Input: nsamples-by-ntimepoints matrix (as outputted by AR.sample.mvn or AR.sample.compositional)
# Returns nothing
############################
AR.plot.samples = function(x #a mtrix as in output of AR.sample.mvn or AR.sample.compositional
){
  y = data.table(x)
  nsamples = nrow(y)
  y[,id := as.factor(1:nsamples)]
  y.m = melt(y, id.vars="id", measure.vars = paste0("V",1:ncol(x)))
  ggplot(y.m, aes(x=as.numeric(variable),y=value, color = id)) + geom_line()
}

