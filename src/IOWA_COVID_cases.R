# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
library( dplyr )
library( tigris )
library( spdep )
library( INLA )
library( tmap )

IA.data = read.csv("./data/IA_COVID_cases.csv")

IA.cum.cases = IA.data[,13:895]   ## days is div by 7
IA.cum.cases.mat = as.matrix( IA.cum.cases )
J = ncol( IA.cum.cases )    ## Number of days before differencing
IA.daily.cases.mat = IA.cum.cases.mat[,2:J] - IA.cum.cases.mat[,1:(J-1)]
J = ncol( IA.daily.cases.mat )
I = nrow( IA.daily.cases.mat )    ## Number of counties
K = J/7                           ## Number of weeks

IA.weekly.cases.mat = matrix( -999 , nrow=I , ncol=K )
for (i in 1:I)
{
  for (k in 1:K)
  {
    IA.weekly.cases.mat[i,k] = sum( IA.daily.cases.mat[i,(7*(k-1)+1):(7*k)] )
  }
}

min( IA.weekly.cases.mat )
for(i in 1:I) print( min( IA.weekly.cases.mat[i,] ) )

IA.weekly.cases.mat.NA = IA.weekly.cases.mat
for ( i in 1:I )
{
  for ( k in 1:K )
  {
    if ( IA.weekly.cases.mat[i,k] < 0 ) IA.weekly.cases.mat.NA[i,k] = NA
  }
}

# windows( 9 , 7 )
plot( 1:K , IA.weekly.cases.mat.NA[1,] , type="l" , ylim=c(-30,3000) , col="gray" )
for (i in 1:I) lines( 1:K , IA.weekly.cases.mat.NA[i,] , col="gray")

all.county.FIPS = as.character( IA.data$countyFIPS )

IA.FIPS = c()
for (i in 1:I)  IA.FIPS = c( IA.FIPS , rep( all.county.FIPS[i],K ) )

IA.cases = c()
for (i in 1:I)  IA.cases = c( IA.cases , IA.weekly.cases.mat.NA[i,] )

IA.county = c()
for (i in 1:I)  IA.county = c( IA.county , rep(i,K) )

IA.week = c()
for (i in 1:I)  IA.week = c( IA.week , 1:K )

IOWA.pop = read.csv( "./data/IOWApop.csv" )$P2020
IA.pop = c()
for (i in 1:I )  IA.pop = c( IA.pop , rep( IOWA.pop[i] , K ) )

IA.df = data.frame( IA.FIPS , IA.county , IA.week , IA.cases , IA.pop )
names( IA.df ) = c( "GEOID" , "county" , "week" , "cases" , "pop" )
IA.df$county1 = IA.county
IA.df$week1   = IA.df$week

#write.csv( IA.df , "IA_df.csv" )


iowa = counties(state = 'IA')   ##  Read in tiger file for geometry
IA.nb = poly2nb( iowa )
nb2INLA( "./dataIOWAmap.adj" , IA.nb )
g = inla.read.graph( filename = "./data/IOWAmap.adj" )


plotIA = function( week )
{
  iowa.week = IA.df[ IA.df$week == week , ]
  P = iowa.week$pop
  y = iowa.week$cases
  COVIDrate = y/P
  iowa.week$COVIDrate = COVIDrate
  iowa.temp = left_join( iowa , iowa.week , by="GEOID" )  
  brks = seq( 0 , 0.0060 , 0.0006 )
  tm_shape( iowa.temp ) + 
    tm_polygons( col = "COVIDrate" , breaks=brks , legend.show = FALSE ) + 
    tm_layout( title = paste0( "Week " , round(week,0) ), 
               title.position = c('left', 'bottom') ,
               title.size = 1 )  
}

p1 = plotIA( 61 )
p2 = plotIA( 62 )
p3 = plotIA( 63 )
p4 = plotIA( 64 )
p5 = plotIA( 65 )
p6 = plotIA( 66 )
p7 = plotIA( 67 )
p8 = plotIA( 68 )
p9 = plotIA( 69 )
p10 = plotIA( 70 )
p11 = plotIA( 71 )
p12 = plotIA( 72 )
p13 = plotIA( 73 )
p14 = plotIA( 74 )
p15 = plotIA( 75 )
p16 = plotIA( 76 )
p17 = plotIA( 77 )
p18 = plotIA( 78 )
p19 = plotIA( 79 )
p20 = plotIA( 80 )
p21 = plotIA( 81 )
p22 = plotIA( 83 )
p23 = plotIA( 84 )
p24 = plotIA( 85 )
p25 = plotIA( 86 )
p26 = plotIA( 87 )
p27 = plotIA( 88 )
p28 = plotIA( 89 )
p29 = plotIA( 90 )
p30 = plotIA( 91 )
p31 = plotIA( 92 )
p32 = plotIA( 93 )
p33 = plotIA( 94 )
p34 = plotIA( 95 )
p35 = plotIA( 96 )

# windows( 15 , 9 )
tmap_arrange( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,
              p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,
              p23,p24,p25,p26,p27,p28,p29,p30,p31,p32,
              p33,p34,p35,
              ncol = 7,nrow = 5 )

formula.par <- cases ~ 1 + f( county ,  model="besag" , graph="./data/IOWAmap.adj" , constr=TRUE ) +
                           f( county1 , model="iid" , constr=TRUE ) +
                           f( week ,    model="rw1" ) +
                           f( week1 ,   model="iid" )

model.par <- inla( formula.par , 
                   family="poisson" , 
                   data=IA.df , 
                   #E=E ,    
                   offset=log(IA.df$pop) ,
                   control.predictor=list(compute=TRUE) ,
                   control.compute=list(dic=TRUE,cpo=TRUE) )

fitted = model.par$summary.fitted.values[,"mean"] 

# windows( 14 , 12 )
par( mfrow=c(4,4) )
i.start = 71
county.i.raw = IA.df$cases[ IA.df$county == i.start ]
county.i.fit = fitted[ IA.df$county == i.start ]
plot(  1:K , county.i.raw , main=paste0("County ",round(i.start)) )
lines( 1:K , county.i.fit )

for ( i in (i.start+1):(i.start+15) )
{
  county.i.raw = IA.df$cases[ IA.df$county == i ]
  county.i.fit = fitted[ IA.df$county == i ]
  plot(  1:K , county.i.raw , main=paste0("County ",round(i)) )
  lines( 1:K , county.i.fit )
}

# windows( 9 , 7 )
i = 71
county.i.raw = IA.df$cases[ IA.df$county == i ]
county.i.fit = fitted[ IA.df$county == i ]
plot(  1:K , county.i.raw , main=paste0("County ",round(i)) )
lines( 1:K , county.i.fit )
text( 30 , 150 , "County 71 has missing counts for weeks 70 & 71" )

################################################################################
####
####  PREDICTION  (Actually, credible interval for future expected value)
####
################################################################################

f = 4   ## number of future weeks to predict

new.week = sort(rep((K+1):(K+f),I))
PREDICT.df = data.frame( GEOID=rep(all.county.FIPS,f) , county=rep(1:I,f) , 
                         week=new.week ,
                         cases=rep(NA,f*I) , pop=rep(IOWA.pop,f) , 
                         county1=rep(1:I,f) , week1=new.week )

IA.df.pred = rbind( IA.df , PREDICT.df ) 

formula.par <- cases ~ 1 + f( county ,  model="besag" , graph="./data/IOWAmap.adj" , constr=TRUE ) +
                           f( county1 , model="iid" , constr=TRUE ) +
                           f( week ,    model="rw1" ) +
                           f( week1 ,   model="iid" )

model.par.pred <- inla( formula.par , 
                        family="poisson" , 
                        data=IA.df.pred , 
                        #E=E ,    
                        offset=log(IA.df.pred$pop) ,
                        control.predictor=list(compute=TRUE,link=1) ,
                        control.compute=list(dic=TRUE,cpo=TRUE) )
#### Note  link=1  in control.predictor above

fitted.pred   = model.par.pred$summary.fitted.values$mean
fitted.pred.L = model.par.pred$summary.fitted.values$`0.025quant`
fitted.pred.U = model.par.pred$summary.fitted.values$`0.975quant`
IA.df.pred$fitted.pred = fitted.pred

# windows( 9 , 7 )
i = 71
county.i.raw = IA.df$cases[ IA.df$county == i ]
county.i.fit = fitted.pred[ IA.df.pred$county == i ]
county.i.fit.L = fitted.pred.L[ IA.df.pred$county == i ]
county.i.fit.U = fitted.pred.U[ IA.df.pred$county == i ]
plot(  1:K , county.i.raw , main=paste0("County ",round(i)) , 
       ylim=c(0,300) , xlim=c(0,K+f) )
lines( 1:(K+f) , county.i.fit )
for ( k in (K+1):(K+f) )
{
  points( k , county.i.fit[k] , col="red" , pch=19 )
  points( k , county.i.fit.L[k] , col="darkred" , pch=2 )
  points( k , county.i.fit.U[k] , col="darkred" , pch=6 )
  lines( c(k,k) , c(county.i.fit.L[k],county.i.fit.U[k]) , col="gray" )
}
points( 70 , county.i.fit[70] , col="red" , pch=19 )
points( 71 , county.i.fit[71] , col="red" , pch=19 )
abline( h=0 )
#text( 30 , 150 , "County 71 has missing counts for weeks 70 & 71" )

################################################################################
####
####  PREDICTION -- Predict a future outcome
####
################################################################################






