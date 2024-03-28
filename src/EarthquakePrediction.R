library( MixtureInf )
library( INLA )
data( earthquake )
earthquake$year = 1900:2006
summary( earthquake )

quake.ar1 = inla( number ~ 1 + f( year , model="ar1" ) ,
                  data=earthquake ,
                  family="poisson" ,
                  control.predictor=list(compute=TRUE) ) 
summary( quake.ar1 )

quake.rw1 = inla( number ~ 1 + f( year , model="rw1" ) ,
                              data=earthquake ,
                              family="poisson" ,
                              control.predictor=list(compute=TRUE) ) 
summary( quake.rw1 )

fitted = quake.ar1$summary.fitted.values
plot( earthquake$year , earthquake$number , pch=19 , col="red" )
lines( earthquake$year , fitted$mean , lwd=2 )
lines( earthquake$year , fitted$`0.025quant` , lwd=1.5 , col="gray" )
lines( earthquake$year , fitted$`0.975quant` , lwd=1.5 , col="gray" )

quake.pred = rbind( earthquake ,
                    data.frame( number=rep(NA,14) , year=2007:2020 ) )

quake.ar1.pred = inla( number ~ 1 + f( year , model="ar1" ) ,
                       data=quake.pred ,
                       family="poisson" ,
                       control.predictor=list(compute=TRUE,link=1) ) 
####  
####  The key to the above is  link=1
####

fitted.pred = quake.ar1.pred$summary.fitted.values
plot( quake.pred$year , quake.pred$number , pch=19 , col="red" , ylim=c(0,50) )
lines( quake.pred$year , fitted.pred$mean , lwd=2 )
lines( quake.pred$year , fitted.pred$`0.025quant` , lwd=1.5 , col="gray" )
lines( quake.pred$year , fitted.pred$`0.975quant` , lwd=1.5 , col="gray" )

