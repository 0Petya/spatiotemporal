
####  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

####  Code taken from
####  https://groups.google.com/g/r-inla-discussion-group/c/tPoKFpK5vqs

##  setwd("C:/Users/srigd/Dropbox/MyCourses/BST_5620_Spring_2023")
# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
library(maptools)
library(sp)
library(INLA)
library(tigris)
library(SpatialEpi)
library(tmap)
library(spdep)
library(nimble)
library(dplyr)

##  url <- "https://sites.google.com/a/r-inla.org/stbook/Ohio_data.zip"
##  temp <- tempfile()
##  download.file(url, temp)
##  unzip(temp, "OhioRespMort1.csv")
OH.data <- read.csv("./data/OhioRespMort.csv")

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

plotLC.Ohio = function( year )
{
  ohio.year = OH.data[ OH.data$year == year , ]
  P = ohio.year$n
  y = ohio.year$y
  LCrate = y/P
  ohio.year$LCrate = LCrate
  ohio.temp = left_join( ohio , ohio.year , by="NAME" )  
  brks = seq( 0 , 0.0014 , 0.0002 )
  tm_shape( ohio.temp ) + 
    tm_polygons( col = "LCrate" , breaks=brks , legend.show = FALSE ) + 
    tm_layout( title = paste0( round(1967+year,0) ), 
               title.position = c('left', 'top') ,
               title.size = 0.8 )  
}

p1 = plotLC.Ohio( 1 )
p2 = plotLC.Ohio( 2 )
p3 = plotLC.Ohio( 3 )
p4 = plotLC.Ohio( 4 )
p5 = plotLC.Ohio( 5 )
p6 = plotLC.Ohio( 6 )
p7 = plotLC.Ohio( 7 )
p8 = plotLC.Ohio( 8 )
p9 = plotLC.Ohio( 9 )
p10 = plotLC.Ohio( 10 )
p11 = plotLC.Ohio( 11 )
p12 = plotLC.Ohio( 12 )
p13 = plotLC.Ohio( 13 )
p14 = plotLC.Ohio( 14 )
p15 = plotLC.Ohio( 15 )
p16 = plotLC.Ohio( 16 )
p17 = plotLC.Ohio( 17 )
p18 = plotLC.Ohio( 18 )
p19 = plotLC.Ohio( 19 )
p20 = plotLC.Ohio( 20 )
p21 = plotLC.Ohio( 21 )

# windows( 15 , 9 )
tmap_arrange( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
              p13,p14,p15,p16,p17,p18,p19,p20,p21,
              ncol = 7,nrow = 3 )

LC.matrix = matrix( 0 , nrow=88 , ncol=21 )    ####  rows are counties; columns are years
for ( year in 1:21 )
{
  ohio.year = OH.data[ OH.data$year == year , ]
  P = ohio.year$n
  y = ohio.year$y
  LCrate = y/P
  LC.matrix[,year] = LCrate
}

# windows( 9 , 7 )
plot( 1:21 , LC.matrix[1,] , col="gray" , ylim=c(0.0001,0.0014) , log="y" )
for (cnty in 2:88 ) lines( 1:21 , LC.matrix[cnty,] , col="gray" )
lines( 1:21 , LC.matrix[18,] , col="red" , lwd=2 )
lines( 1:21 , LC.matrix[1,] , col="blue" , lwd=2 )
lines( 1:21 , LC.matrix[57,] , col="green" , lwd=2 )

OH.data$LCrate = OH.data$y / OH.data$n

################################################################################
####
####  MCMC for no-interaction model
####
################################################################################

I = max( OH.data$county )    ## for OH data, I = 88 counties
J = max( OH.data$year )      ## for OH data, J = 21 time periods (years)
y.matrix = matrix( 0 , nrow=I , ncol=J )
P.matrix = matrix( 0 , nrow=I , ncol=J )

for ( j in 1:J )
{
  ohio.year = OH.data[ OH.data$year == j , ]
  y.matrix[,j] = ohio.year$y
  P.matrix[,j] = ohio.year$n
}

spatialCode = nimbleCode({
  
  ##  Likelihood
  ##  
  for (i in 1:I) 
  {
    for (j in 1:J)
    {
      y[i,j] ~ dpois( P[i,j]*theta[i,j] )
      log(theta[i,j])  <-  beta0 + beta1*t[j] + s[i] + v[i]
    }
  }
  
  ##  Latent Variables
  for (l in 1:L)
  {
    weights[l]  <-  1
  }
  s[1:I] ~ dcar_normal(adj[1:L],weights[1:L],num[1:I],taus,zero_mean=1)
  for (i in 1:I)
  {
    v[i] ~ dnorm(0,tauv)
  }
  
  ##  Hyperpriors
  beta0 ~ dnorm( 0 , 0.01 )
  beta1 ~ dnorm( 0 , 0.01 )
  taus ~ dgamma( 1 , 0.001 )
  tauv ~ dgamma( 1 , 0.001 )
})

constants = list( I=I , J=J , L=L , num=num , adj=adj , P=P.matrix , t=(1:J) )
data = list( y=y.matrix )
inits = list( beta0=-8 , beta1=0 , taus=1 , tauv = 1 , s=rep(0,I) , v=rep(0,I) )

spatial.model = nimbleModel( code=spatialCode , 
                             constants=constants , 
                             data=data ,
                             inits=inits )

compile.spatial.model = compileNimble( spatial.model )

spatial.model.Conf = configureMCMC( spatial.model , print = TRUE , thin=100 )

spatial.model.Conf$addMonitors( c( "beta0" , "beta1" , "taus" , "tauv" , "theta" ) )

spatial.model.MCMC = buildMCMC( spatial.model.Conf )

compile.spatial.MCMC = compileNimble( spatial.model.MCMC )

niter =  50000
nburn =  20000
#set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

str( samples.spatial )

beta0.samp = as.vector(samples.spatial[,"beta0"])
beta1.samp = as.vector(samples.spatial[,"beta1"])
taus.samp = as.vector(samples.spatial[,"taus"])
tauv.samp = as.vector(samples.spatial[,"tauv"])
theta.1.1.samp = as.vector(samples.spatial[,"theta[1, 1]"])
theta.1.21.samp = as.vector(samples.spatial[,"theta[1, 21]"])

# windows( 9 , 12 )
par( mfrow=c(3,2) )
ts.plot( beta0.samp , main="Trace Plot for beta0" )
ts.plot( beta1.samp , main="Trace Plot for beta1" )
ts.plot( taus.samp , main="Trace Plot for taus" )
ts.plot( tauv.samp , main="Trace Plot for tauv" )
ts.plot( theta.1.1.samp , main="Trace Plot for theta[1,1]" )
ts.plot( theta.1.21.samp , main="Trace Plot for theta[1,21]" )

# windows( 6 , 6 )
plot( taus.samp , tauv.samp )

head( samples.spatial )
theta.hat    = apply( samples.spatial[,5:((I*J)+4)] , 2 , mean )
se.theta.hat = apply( samples.spatial[,5:((I*J)+4)] , 2 , sd )

LCL = theta.hat - 1.96*se.theta.hat
UCL = theta.hat + 1.96*se.theta.hat

# windows( 14 , 10 )
par( mfrow=c(4,6))
for( county in 49:72)
{
  plot( 1:J , log(LC.matrix[county,]) , ylim=c(-10,-6) , 
        ylab="log(theta)" , xlab="Year" ,
        main=paste0(OH.data$NAME[county]," County ") )
  lines( 1:J , log( theta.hat[seq(from=county,by=I,length=J)] ) )
}

# windows( 14 , 10 )
par( mfrow=c(4,6))
for( county in 49:72)
{
  plot( 1:J , LC.matrix[county,] , ylim=c(0,0.0014) , 
        ylab="theta" , xlab="Year" ,
        main=paste0(OH.data$NAME[county]," County ") )
  lines( 1:J , theta.hat[seq(from=county,by=I,length=J)] )
}




