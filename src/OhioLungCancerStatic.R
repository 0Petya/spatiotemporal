
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
##  unzip(temp, "Ohio_data/OhioRespMort.csv")
OH.data <- read.csv("OhioRespMort.csv")

##  unzip(temp, exdir = tempdir())
##  shapefile <- paste0(tempdir(),"\\Ohio_data\\tl_2010_39_county00")
##  ohio <- readShapeSpatial(shapefile)
##  ohio = rgdal::readOGR(tl_2010_39_county00)

ohio <- counties(state = 'OH')    ##  Get shape file from tigris
View( ohio )                      ##  Let's investigate the object

class( ohio )
ohio.nb = poly2nb( ohio )
ohio.net = nb2lines( nc.sids.nb , coords=coordinates(nc.sids) )

windows( 14 , 8 )
tm_shape( ohio ) +
  tm_borders( col="darkgray" ) +
  tm_shape( ohio.net ) +
  tm_lines( col="darkgreen" , lwd=2 )

ohio.lw = nb2listw( ohio.nb )

k = nrow( ohio )
num = rep(0,k)
for (i in 1:k) num[i] = length( ohio.lw$neighbours[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,ohio.lw$neighbours[[i]] )
L = length(adj)

ohio.year1 = OH.data[ OH.data$year == 1 , ]

names( ohio.year1 )
y = ohio.year1$y
P = ohio.year1$n
E = ohio.year1$E
R = sum( y ) / sum( P )
E1 = P*R                        ## Note E from file agrees with E1 computed here

spatialCode = nimbleCode({
  mu ~ dflat()
  tau ~ dgamma( 1 , 0.001 )
  for (i in 1:L)
    weights[i]  <-  1
  s[1:k] ~ dcar_normal(adj[1:L],weights[1:L],num[1:k],tau,zero_mean=1)
  for (i in 1:k) {
    log(theta[i])  <-  mu + s[i]
    y[i] ~ dpois( E[i]*theta[i] )
  }  
})

constants = list( k=k , L=L , num=num , adj=adj , E=E )
data = list( y=y )
inits = list( mu=0 , tau=1 , s=rep(0,k) )

ohio.spatial.model = nimbleModel( code=spatialCode , 
                                  constants=constants , 
                                  data=data ,
                                  inits=inits )

compile.ohio.spatial.model = compileNimble( ohio.spatial.model )

ohio.spatial.model.Conf = configureMCMC( ohio.spatial.model , print = TRUE )

ohio.spatial.model.Conf$addMonitors(c("mu","tau","theta"))

ohio.spatial.model.MCMC = buildMCMC( ohio.spatial.model.Conf )

compile.ohio.spatial.MCMC = compileNimble( ohio.spatial.model.MCMC )

niter = 21000
nburn =  1000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.ohio.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

str( samples.spatial )
head( samples.spatial )

theta.hat    = apply( samples.spatial[,3:90] , 2 , mean )
se.theta.hat = apply( samples.spatial[,3:90] , 2 , sd )
LCL = theta.hat - 1.96*se.theta.hat
UCL = theta.hat + 1.96*se.theta.hat

windows( 14 , 7 )
plot( 1:88 , theta.hat , ylim=c(0,2) , xlab("County") , ylab="SMR" )
for (i in 1:88 ) lines( c(i,i) , c(LCL[i],UCL[i]) )
abline( h = 1 )
