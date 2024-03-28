# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
GA.data = read.csv( "./data/GeorgiaEmergencyRoomVIsits.csv" )[1:159,]
GA.data$NAME = GA.data$County
GA.pop = read.csv( "./data/GeorgiaPopulation.csv" , header=FALSE )
names( GA.pop ) = c( "NAME" , "pop" )

year = 2002:2021
cnty = 56
# windows( 9 , 7 )
plot( year , as.vector(GA.data[cnty,2:21]) , type="b" , ylab="Emer Room Visits" )

library(maptools)
library(sp)
library(INLA)
library(tigris)
library(SpatialEpi)
library(tmap)
library(spdep)
library(nimble)
library(dplyr)

georgia <- counties(state = 'GA')    ##  Get shape file from tigris

georgia.nb = poly2nb( georgia )
georgia.lw = nb2listw( georgia.nb )

k = nrow( georgia )
num = rep(0,k)
for (i in 1:k) num[i] = length( georgia.lw$neighbours[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,georgia.lw$neighbours[[i]] )
L = length(adj)


GA.data = left_join( GA.data , GA.pop )

plotERV.georgia = function( year )
{
  colmn = paste0("Y",round(year))
  GA.temp = GA.data
  P = GA.data$pop
  y = GA.data[,colmn]
  ERVrate = y/P
  GA.temp$ERVrate = ERVrate
  georgia.temp = left_join( georgia , GA.temp , by="NAME" )  
  brks = seq( 0 , 1 , 0.1 )
  tm_shape( georgia.temp ) + 
    tm_polygons( col = "ERVrate" , breaks=brks , legend.show = FALSE ) + 
    tm_layout( title = paste0( round(year,0) ), 
               title.position = c('right', 'top') ,
               title.size = 1 )  
}

p1 = plotERV.georgia( 2002 )
p2 = plotERV.georgia( 2003 )
p3 = plotERV.georgia( 2004 )
p4 = plotERV.georgia( 2005 )
p5 = plotERV.georgia( 2006 )
p6 = plotERV.georgia( 2007 )
p7 = plotERV.georgia( 2008 )
p8 = plotERV.georgia( 2009 )
p9 = plotERV.georgia( 2010 )
p10 = plotERV.georgia( 2011 )
p11 = plotERV.georgia( 2012 )
p12 = plotERV.georgia( 2013 )
p13 = plotERV.georgia( 2014 )
p14 = plotERV.georgia( 2015 )
p15 = plotERV.georgia( 2016 )
p16 = plotERV.georgia( 2017 )
p17 = plotERV.georgia( 2018 )
p18 = plotERV.georgia( 2019 )
p19 = plotERV.georgia( 2020 )
p20 = plotERV.georgia( 2021 )

# windows( 14 , 10 )
tmap_arrange( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
              p13,p14,p15,p16,p17,p18,p19,p20,
              ncol = 5,nrow = 4 )

ERV.matrix = matrix( 0 , nrow=159 , ncol=20 )    ####  rows are counties; columns are years

for ( i in 1:20 )
{
  colmn = paste0( "Y" , round(i+2001) )
  ERV.matrix[,i] = GA.data[,colmn] / GA.data[,"pop"]
}

# windows( 9 , 7 )
plot( 1:20 , ERV.matrix[1,] , col="gray" , type="l" , ylim=c(0,1) )
for (cnty in 2:88 ) lines( 1:20 , ERV.matrix[cnty,] , col="gray" )
lines( 1:20 , ERV.matrix[16,] , col="red" , lwd=2 )      ##  Bulloch
lines( 1:20 , ERV.matrix[44,] , col="blue" , lwd=2 )     ##  DeKalb
lines( 1:20 , ERV.matrix[60,] , col="purple" , lwd=2 )   ##  Fulton
lines( 1:20 , ERV.matrix[131,] , col="green" , lwd=2 )   ##  Taliaferro

################################################################################
####
####  MCMC for First-Order Random Walk
####
################################################################################


I = 159                      ## for GA data, I = 159 counties
J = 20                       ## for GA ERV data, J = 20 time periods (years)
y.matrix = data.matrix( GA.data[,2:21] )
P.vector = GA.pop$pop

spatialCode = nimbleCode({
  
  ##  Likelihood
  ##  
  for (i in 1:I) 
  {
    for (j in 1:J)
    {
      y[i,j] ~ dpois( P[i]*theta[i,j] )
      log(theta[i,j])  <-  beta0 + gamma[j] + phi[j] + s[i] + v[i]
    }
  }
  
  ##  Latent Variables
  gamma[1] ~ dnorm( 0 , 1 )               ##  Initial gamma must be given a prior
  for (j in 2:J)                          ##  Note: loop starts at 2
  {
    gamma[j] ~ dnorm( gamma[j-1] , taug )
  }
  for (j in 1:J)
  {
    phi[j] ~ dnorm( 0 , taup )
  }
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
  taug ~ dgamma( 1 , 0.001 )
  taup ~ dgamma( 1 , 0.001 )
  taus ~ dgamma( 1 , 0.001 )
  tauv ~ dgamma( 1 , 0.001 )
})


constants = list( I=I , J=J , L=L , num=num , adj=adj , P=P.vector , t=(1:J) )
data = list( y=y.matrix )
inits = list( beta0=-2 , taug=1 , taup = 1 , taus=1 , tauv = 1 , 
              s=rep(0,I) , v=rep(0,I) )

spatial.model = nimbleModel( code=spatialCode , 
                             constants=constants , 
                             data=data ,
                             inits=inits )

compile.spatial.model = compileNimble( spatial.model )

spatial.model.Conf = configureMCMC( spatial.model , print = TRUE , thin=100 )

spatial.model.Conf$addMonitors( c( "beta0" , "gamma" , "phi" , "taug" , "taup" , "taus" , 
                                   "tauv" , "theta" ) )

spatial.model.MCMC = buildMCMC( spatial.model.Conf )

compile.spatial.MCMC = compileNimble( spatial.model.MCMC )

niter =  25000
nburn =   5000
#set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )



head( samples.spatial )

beta0.samp = as.vector(samples.spatial[,"beta0"])
taug.samp = as.vector(samples.spatial[,"taug"])
taup.samp = as.vector(samples.spatial[,"taup"])
taus.samp = as.vector(samples.spatial[,"taus"])
tauv.samp = as.vector(samples.spatial[,"tauv"])
gamma.10.samp = as.vector(samples.spatial[,"gamma[10]"])
phi.10.samp   = as.vector(samples.spatial[,"phi[10]"])
theta.1.1.samp = as.vector(samples.spatial[,"theta[1, 1]"])
theta.1.20.samp = as.vector(samples.spatial[,"theta[1, 20]"])

# windows( 9 , 8 )
plot( gamma.10.samp , phi.10.samp )

# windows( 9 , 12 )
par( mfrow=c(4,2) )
ts.plot( beta0.samp , main="Trace Plot for beta0" )
ts.plot( taus.samp , main="Trace Plot for taus" )
ts.plot( tauv.samp , main="Trace Plot for tauv" )
ts.plot( taud.samp , main="Trace Plot for taud" )
ts.plot( gamma.10.samp , main="Trace Plot for gamma[10]" )
ts.plot( phi.10.samp , main="Trace Plot for phi[10]" )
ts.plot( theta.1.1.samp , main="Trace Plot for theta[1,1]" )
ts.plot( theta.1.20.samp , main="Trace Plot for theta[1,20]" )

theta.hat = apply( samples.spatial[,46:3225] , 2 , mean )
gamma.hat = apply( samples.spatial[,2:21] , 2 , mean )

theta.hat.matrix = matrix( theta.hat , nrow=I , ncol=J )

# windows( 12 , 9 )
par( mfrow=c(4,6) , mar=c(5,4,3,2) )
for (i in 41:64)
{
  plot( 1:J , ERV.matrix[i,] , ylim=c(0,1) , main=GA.data$County[i] )
  lines( 1:J , theta.hat.matrix[i,] )
}

