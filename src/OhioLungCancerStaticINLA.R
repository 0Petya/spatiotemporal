
####  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

####  Code taken from
####  https://groups.google.com/g/r-inla-discussion-group/c/tPoKFpK5vqs

setwd("C:/Users/srigd/Dropbox/MyCourses/BST_5620_Spring_2023")
setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
library(maptools)
library(sp)
library(INLA)
library(tigris)
library(SpatialEpi)
library(tmap)
library(spdep)
library(nimble)

##  url <- "https://sites.google.com/a/r-inla.org/stbook/Ohio_data.zip"
##  temp <- tempfile()
##  download.file(url, temp)
##  unzip(temp, "OhioRespMort1.csv")
OH.data <- read.csv("OhioRespMort.csv")

##  unzip(temp, exdir = tempdir())
##  shapefile <- paste0(tempdir(),"\\Ohio_data\\tl_2010_39_county00")
##  ohio <- readShapeSpatial(shapefile)
##  ohio = rgdal::readOGR(tl_2010_39_county00)

ohio <- counties(state = 'OH')    ##  Get shape file from tigris
View( ohio )                      ##  Let's investigate the object

class( ohio )
ohio.nb = poly2nb( ohio )

ohio.lw = nb2listw( ohio.nb )

k = nrow( ohio )
num = rep(0,k)
for (i in 1:k) num[i] = length( ohio.lw$neighbours[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,ohio.lw$neighbours[[i]] )
L = length(adj)

ohio.year1 = OH.data[ OH.data$year == 1 , ]

names( ohio.year1 )
P = ohio.year1$n
y = ohio.year1$y

spatialCode = nimbleCode({
  for (i in 1:k) {
    y[i] ~ dpois( P[i]*theta[i] )
  }  
  for (i in 1:k) {
    log(theta[i])  <-  mu + s[i] + v[i]
  }    
  for (i in 1:L)
  {
    weights[i]  <-  1
  }
  s[1:k] ~ dcar_normal( adj[1:L] , weights[1:L] , num[1:k] , taus , zero_mean=1 )
  for (i in 1:k)
  {
    v[i] ~ dnorm( 0 , tauv )
  }
  mu ~ dnorm(0,0.01)
  taus ~ dgamma( 1 , 0.01 )
  tauv ~ dgamma( 1 , 0.01 )
})

R = sum( ohio.year1$y ) / sum( ohio.year1$n )

constants = list( k=k , L=L , num=num , adj=adj , P=P )
data = list( y=y )
inits = list( mu=0 , taus=1 , tauv=1 , s=rep(0,k) , v=rep(0,k) , theta=rep(R,k) )

spatial.model = nimbleModel( code=spatialCode , 
                                  constants=constants , 
                                  data=data ,
                                  inits=inits )

compile.spatial.model = compileNimble( spatial.model )

spatial.model.Conf = configureMCMC( spatial.model , print = TRUE )

spatial.model.Conf$addMonitors(c("mu","taus","tauv","theta"))

spatial.model.MCMC = buildMCMC( spatial.model.Conf )

compile.spatial.MCMC = compileNimble( spatial.model.MCMC )

niter = 55000
nburn =  5000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

str( samples.spatial )
head( samples.spatial )

windows( 9 , 6 )
plot( samples.spatial[,4] )

tauv.sample = as.vector( samples.spatial[,"tauv"] )
taus.sample = as.vector( samples.spatial[,"taus"] )
plot( tauv.sample , taus.sample , log="xy" )

tauvhat = mean( tauv.sample )
taushat = mean( taus.sample )

#############################################################
####
####  Using INLA
####
#############################################################

ohio.nb = poly2nb( ohio )
nb2INLA( "ohio.adj" , ohio.nb )
g = inla.read.graph(filename="ohio.adj")    ##  Read the graph right back in

ohio.year1$idarea  = 1:k                    ##  INLA needs different IDs for county, one
ohio.year1$idarea1 = 1:k                    ##    for the  besag  and one for the  iid.

hyperPrior1 = list( prec = list(prior="loggamma",param=c(1,0.01)) )
hyperPrior2 = list( prec = list(prior="loggamma",param=c(1,0.01)) )

formula1 =  y ~ f( idarea , model="besag" , graph=g , hyper=hyperPrior1 ) +
                f( idarea1 , model="iid"  , graph=g , hyper=hyperPrior2 )
mod1 = inla( formula1 , family="poisson" , data=ohio.year1 , offset=log(n) )
summary( mod1 )
View( mod1 )

formula2 =  y ~ f( idarea , model="bym" , graph=g , hyper=hyperPrior1 ) 
mod2 = inla( formula2 , family="poisson" , data=ohio.year1 , offset=log(n) )
summary( mod2 )
View( mod2 )

#############################################################
####
####  BYM2 Using INLA
####
#############################################################

formula3 =  y ~ f( idarea , model="bym2" , graph=g ) 
mod3 = inla( formula3 , family="poisson" , data=ohio.year1 , offset=log(n) )
summary( mod3 )
View( mod3 )

#############################################################
####
####  BYM2 Using MCMC and Nimble
####
#############################################################

spatialCode = nimbleCode({
  for (i in 1:k) {
    y[i] ~ dpois( P[i]*theta[i] )
  }  
  for (i in 1:k) {
    log(theta[i])  <-  mu + b[i]
  }    
  for (i in 1:k) {
    b[i] <- ( sqrt(1-phi)*v[i] + sqrt(phi)*s[i] ) / sqrt(taub)
  } 
  for (i in 1:L)
  {
    weights[i]  <-  1
  }
  s[1:k] ~ dcar_normal( adj[1:L] , weights[1:L] , num[1:k] , 1 , zero_mean=1 )
  for (i in 1:k)
  {
    v[i] ~ dnorm( 0 , 1 )
  }
  mu ~ dnorm(0,0.01)
  taub ~ dgamma( 1 , 0.01 )
  phi  ~ dunif(0,1)
})

R = sum( ohio.year1$y ) / sum( ohio.year1$n )

constants = list( k=k , L=L , num=num , adj=adj , P=P )
data = list( y=y )
inits = list( mu=0 , taub=1 , phi=0.5 , s=rep(0,k) , v=rep(0,k) , theta=rep(R,k) )

spatial.model = nimbleModel( code=spatialCode , 
                             constants=constants , 
                             data=data ,
                             inits=inits )

compile.spatial.model = compileNimble( spatial.model )

spatial.model.Conf = configureMCMC( spatial.model , print = TRUE )

spatial.model.Conf$addMonitors(c("mu","taub","phi","theta"))

spatial.model.MCMC = buildMCMC( spatial.model.Conf )

compile.spatial.MCMC = compileNimble( spatial.model.MCMC )

niter = 55000
nburn =  5000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

str( samples.spatial )
head( samples.spatial )

hist( samples.spatial[,"phi"] )
