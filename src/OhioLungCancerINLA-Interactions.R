
####  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

####  Code taken from
####  https://groups.google.com/g/r-inla-discussion-group/c/tPoKFpK5vqs
# setwd("C:/Users/srigd/Dropbox/OhioLungCancer")
# setwd("C:/Users/serig/Dropbox/OhioLungCancer")
# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5620_Spring_2023")
library(maptools)
library(sp)
library(INLA)
library(tigris)
library(SpatialEpi)
library(tmap)
library(dplyr)

##  url <- "https://sites.google.com/a/r-inla.org/stbook/Ohio_data.zip"
##  temp <- tempfile()
##  download.file(url, temp)
##  unzip(temp, "Ohio_data/OhioRespMort.csv")
OH.data <- read.csv("./data/OhioRespMort.csv")

##  unzip(temp, exdir = tempdir())
##  shapefile <- paste0(tempdir(),"\\Ohio_data\\tl_2010_39_county00")
##  ohio <- readShapeSpatial(shapefile)
##  ohio = rgdal::readOGR(tl_2010_39_county00)

ohio <- counties(state = 'OH')
county1 <- OH.data$county
Ohio.adj <- inla.read.graph("./data/ohio.graph")

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

################################################################################
###
###   BYM Model:  s & v   along with linear trend
###
################################################################################
OH.data$county1 = OH.data$county 

formula.par <- y ~ 1 + f( county ,  model="besag" , graph=Ohio.adj, constr=TRUE ) +
                       f( county1 , model="iid" , constr=TRUE ) +
                       year

model.par <- inla( formula.par , 
                   family="poisson" , 
                   data=OH.data , 
                   #E=E ,    
                   offset=log(OH.data$n) ,
                   control.predictor=list(compute=TRUE) ,
                   control.compute=list(dic=TRUE,cpo=TRUE) )

fitted = model.par$summary.fitted.values[,"mean"] 

model.par$summary.fixed[,1:5]

x  = model.par$marginals.random
x1 = model.par$marginals.random[[1]]
x2 = model.par$marginals.random[[2]]

# windows( 9 , 9 )
par( mfrow=c(2,1) )
plot( x1$index.1[,1] , x1$index.1[,2] , type="b" )
plot( x2$index.1[,1] , x2$index.1[,2] , type="b" )

I = 88
J = 21
str( fitted )
OH.LC.fitted.Linear = matrix( fitted , nrow=I , ncol=J )
Pop = matrix( OH.data$n , nrow=I , ncol=J )
OH.LC.fitted.Rate.Linear = OH.LC.fitted.Linear / Pop

# windows( 9 , 7 )
plot( 1:J , OH.LC.fitted.Rate.Linear[1,] , type="l" , col="gray" , ylim=c(0,0.0008) )
for( i in 2:I ) lines( 1:J , OH.LC.fitted.Rate.Linear[i,] , col="gray" )
lines( 1:J , OH.LC.fitted.Rate.Linear[57,] , col="red" , lwd=2 )

################################################################################
###
###   BYM Model:  s & v   along with 1st-order random walk
###
################################################################################

OH.data$year1 = OH.data$year

formula.par <- y ~ 1 + f( county ,  model="besag" , graph=Ohio.adj, constr=TRUE ) +
                       f( county1 , model="iid" , constr=TRUE ) +
                       f( year ,    model="rw1" ) +
                       f( year1 ,   model="iid" )

model.par <- inla( formula.par , 
                   family="poisson" , 
                   data=OH.data , 
                   #E=E ,    
                   offset=log(OH.data$n) ,
                   control.predictor=list(compute=TRUE) ,
                   control.compute=list(dic=TRUE,cpo=TRUE) )

fitted = model.par$summary.fitted.values[,"mean"] 

model.par$summary.fixed[,1:5]

x  = model.par$marginals.random
x1 = model.par$marginals.random[[1]]
x2 = model.par$marginals.random[[2]]

# windows( 9 , 9 )
par( mfrow=c(2,1) )
plot( x1$index.1[,1] , x1$index.1[,2] , type="b" )
plot( x2$index.1[,1] , x2$index.1[,2] , type="b" )

I = 88
J = 21
str( fitted )
OH.LC.fitted.RW1 = matrix( fitted , nrow=I , ncol=J )
Pop = matrix( OH.data$n , nrow=I , ncol=J )
OH.LC.fitted.Rate.RW1 = OH.LC.fitted.RW1 / Pop

# windows( 9 , 7 )
plot( 1:J , OH.LC.fitted.Rate.RW1[1,] , type="l" , col="gray" , ylim=c(0,0.0008) )
for( i in 2:I ) lines( 1:J , OH.LC.fitted.Rate.RW1[i,] , col="gray" )
lines( 1:J , OH.LC.fitted.Rate.RW1[57,] , col="red" , lwd=2 )

################################################################################
###
###   BYM Model:  s & v   along with 1st-order random walk and 
###               Type I interaction (phi & v)
###
################################################################################

OH.data$area.year = 1:(I*J)

formula.par <- y ~ 1 + f( county ,  model="besag" , graph=Ohio.adj, constr=TRUE ) +
                       f( county1 , model="iid" , constr=TRUE ) +
                       f( year ,    model="rw1" ) +
                       f( year1 ,   model="iid" ) +
                       f( area.year , model="iid" )
  
model.I  <-  inla( formula.par , 
                   family="poisson" , 
                   data=OH.data , 
                   #E=E ,    
                   offset=log(OH.data$n) ,
                   control.predictor=list(compute=TRUE) ,
                   control.compute=list(dic=TRUE,cpo=TRUE) )

fitted = model.I$summary.fitted.values[,"mean"] 

model.I$summary.fixed[,1:5]

I = 88
J = 21
str( fitted )
OH.LC.fitted.TypeI = matrix( fitted , nrow=I , ncol=J )
Pop = matrix( OH.data$n , nrow=I , ncol=J )
OH.LC.fitted.Rate.TypeI = OH.LC.fitted.TypeI / Pop

# windows( 9 , 7 )
plot( 1:J , OH.LC.fitted.Rate.TypeI[1,] , type="l" , col="gray" , 
      main="Type I Interaction" , ylim=c(0,0.0008) )
for( i in 2:I ) lines( 1:J , OH.LC.fitted.Rate.TypeI[i,] , col="gray" )
lines( 1:J , OH.LC.fitted.Rate.TypeI[57,] , col="red" , lwd=2 )

################################################################################
###
###   BYM Model:  s & v   along with 1st-order random walk and 
###               Type II interaction (phi & v)
###
################################################################################

OH.data$area.int = OH.data$county
OH.data$year.int = OH.data$year
formula.II =   y  ~  f( county , model="besag" , graph=Ohio.adj, constr=TRUE ) +
                     f( county1 , model="iid" , constr=TRUE ) +
                     f( year ,    model="rw1" ) +
                     f( year1 ,   model="iid" ) +
                     f( area.int , model="iid" , group=year.int ,
                        control.group=list(model="rw1") )
model.II <- inla( formula.II , 
                   family="poisson" , 
                   data=OH.data , 
                   #E=E ,    
                   offset=log(OH.data$n) ,
                   control.predictor=list(compute=TRUE) ,
                   control.compute=list(dic=TRUE,cpo=TRUE) )

fitted.II = model.II$summary.fitted.values[,"mean"] 

model.I$summary.fixed[,1:5]

model.I$summary.fixed[,1:5]

I = 88
J = 21
str( fitted.II )
OH.LC.fitted.TypeII = matrix( fitted , nrow=I , ncol=J )
Pop = matrix( OH.data$n , nrow=I , ncol=J )
OH.LC.fitted.Rate.TypeII = OH.LC.fitted.TypeII / Pop

# windows( 9 , 7 )
plot( 1:J , OH.LC.fitted.Rate.TypeII[1,] , type="l" , col="gray" , 
                       main="Type II Interaction" , ylim=c(0,0.0008) )
for( i in 2:I ) lines( 1:J , OH.LC.fitted.Rate.TypeII[i,] , col="gray" )
lines( 1:J , OH.LC.fitted.Rate.TypeII[57,] , col="red" , lwd=2 )
